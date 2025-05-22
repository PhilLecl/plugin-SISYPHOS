# The .sisy jobs file contains all the information neccecary to perform a iterative refinement of a given structure.
# Defined are mostly parameters that are later supplied to NoSpherA2 via the olex2c-headless or olex2-gui version
# The .sisy file follows these rules:
# Every line containing # or only /n (a new line) is ignored
# Every other line should contain a single job as follows:
# NoSpherA2_Keyword1:Value1;Keyword2:Value2;Keyword3:Value3... and so on.
# When reading each line we split the line at the ; and then at the first(!): to get a dictionary of keyword: value pairs (We split at the first : to allow for keywords containing : e.g. File paths under windows)
import glob
import itertools as it
import os
from typing import Dict, List, Tuple


class SisyphosBenchmarkFile:
    def __init__(self, location: str, options=None):
        self.jobfile = glob.glob(os.path.join(location, "*.sisy"))[0]
        self.location = location
        self.jobs: list[dict[str, str]] = []
        if options == None:
            self.read()
        else:
            self.write(options)

    def read(self):
        with open(self.jobfile, "r") as f:
            for line in f:
                if line.startswith("#") or line == "\n":
                    continue
                job = {}
                for item in line.split(";"):
                    if ":" in item:
                        key, value = item.split(":", maxsplit=1)
                        job[key.strip()] = value.strip()
                self.jobs.append(job)

    def __len__(self):
        return len(self.jobs)

    def __getitem__(self, idx):
        return self.jobs[idx]

    def __iter__(self):
        for job in self.jobs:
            yield job

    def is_finished(self, job_id: int) -> bool:
        """Check if the job with the given id is finished"""
        jobdir = os.path.join(self.location, f"job_{job_id}")
        done_file = os.path.join(jobdir, "done")
        if not os.path.exists(done_file):
            return False
        return True

    def check_all_done(self) -> Tuple[List[int], List[int]]:
        """Iterate through all directorys and check if the 'done' file is present
        Returns a tuple of two lists, the first one contains the finished jobs, the second one the running jobs
        """
        finished = []
        running = []
        for job_id in range(len(self.jobs)):
            if self.is_finished(job_id):
                finished.append(job_id)
            else:
                running.append(job_id)

        return (finished, running)

    def write(self, OPTIONS: Dict[str, str]) -> int:
        permutations = [{"IAM": "True"}]
        # As SALTED does not rely on ORCA calculation, throw out everything that is not needed
        if OPTIONS["selected_salted_model"]:
            for model in OPTIONS["selected_salted_model"]:
                permutations.append(
                    {"selected_salted_model": model, "source": "SALTED"}
                )

        # Now calculte the permutations for the rest of the options
        del OPTIONS["selected_salted_model"]
        for perm in it.product(*OPTIONS.values()):
            temp_dict = dict(zip(OPTIONS.keys(), perm))
            if temp_dict["NoSpherA2_RI_FIT"] == "False":
                # If NoSpherA2_RI_FIT is False, we do not need to permutate the auxiliary_basis and auto_aux_beta as it is not used
                if (temp_dict["auxiliary_basis"] != OPTIONS["auxiliary_basis"][0]) or (
                    temp_dict["auto_aux_beta"] != OPTIONS["auto_aux_beta"][0]
                ):
                    continue
            else:
                # If NoSpherA2_RI_FIT is True, but auto_aux is not used, we do not need to permutate the beta value
                if (temp_dict["auxiliary_basis"] != "auto_aux") and (
                    temp_dict["auto_aux_beta"] != OPTIONS["auto_aux_beta"][0]
                ):
                    continue

            permutations.append(temp_dict)

        # Write the permutations to the file
        with open(os.path.join(self.jobfile), "w") as f:
            f.write(f"#SISYPHOS creating {len(permutations)} tasks!\n")
            f.write("#The following options are available:\n")
            for key in OPTIONS.keys():
                f.write(f"# {key}: {OPTIONS[key]}\n")
            f.write("#\n")
            f.write("#The following options are used:\n")
            for perm in permutations:
                n_keys = len(perm)
                for i, (key, value) in enumerate(perm.items()):
                    f.write(f"{key}:{value}")
                    if i < n_keys - 1:
                        f.write(";")
                f.write("\n")

        return len(permutations)
