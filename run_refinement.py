# PLACE THIS FILE IN THE OLEX2 FOLDER!
import glob
import itertools as it
import multiprocessing
import os
import subprocess
import sys
from typing import Tuple

# GLOBAL CONFIGURATION OF SISIPHOS!
ncpus = 16
cpus_per_job = 4


OPTIONS = {
    "basis_name": ["def2-SVP", "def2-TZVP"],
    "method": ["r2SCAN"],
    "charge": ["0"],
    "multiplicity": ["0"],
    "becke_accuracy": ["Normal"],
    "h_aniso": ["True"],
    "h_afix": ["False"],
    "ORCA_SCF_Conv": ["NoSpherA2SCF"],
    "ORCA_SCF_Strategy": ["EasyConv"],
    "ORCA_Solvation": ["Water"],
    "NoSpherA2_RI_FIT": ["True", "False"],
    "auxiliary_basis": ["auto_aux"],
    "auto_aux_beta": ["1.5", "1.7", "2.0", "2.5"],
    "selected_salted_model": [
        R"D:\Models\Combo_v1",
        R"D:\Models\Combo_v2",
        R"D:\Models\Combo_v3",
        R"D:\Models\Combo_v4",
        R"D:\Models\Combo_v5",
        R"D:\Models\Combo_v6",
    ],
}


def write_sisy_file(work_dir: str, filename: str) -> int:
    permutations = [{"IAM": "True"}]
    # As SALTED does not rely on ORCA calculation, throw out everything that is not needed
    if OPTIONS["selected_salted_model"]:
        for model in OPTIONS["selected_salted_model"]:
            permutations.append({"selected_salted_model": model, "source": "SALTED"})

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
    with open(os.path.join(work_dir, filename), "w") as f:
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


def execute_olex2c(olex_exe: str, jobdir: str, job_idx: int) -> None:
    my_env = os.environ.copy()
    my_env["PYTHONHOME"] = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "Python38"
    )
    my_env["PYTHONPATH"] = os.path.join(my_env["PYTHONHOME"], "Lib")
    my_env["SISYPHOS_job_idx"] = str(job_idx)

    with open(os.path.join(jobdir, f"output_{job_idx}.log"), "w") as f:
        f.write(os.path.dirname(os.path.abspath(__file__)))
        try:
            process = subprocess.Popen(
                args=[olex_exe, "sisyphos"],
                stdin=subprocess.PIPE,
                stdout=f,
                stderr=f,
                env=my_env,
                cwd=os.path.dirname(olex_exe),
            )
            process.communicate()
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, "olex2c.exe")
        except subprocess.CalledProcessError as e:
            f.write(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
        except Exception as e:
            f.write(f"Error in execute_olex2c: {e}")

    # The outfile is written in different encodings, this is a workaround to fix it
    with open(os.path.join(jobdir, f"output_{job_idx}.log"), "rb") as f:
        data = f.read().replace(b"\x00", b"").decode("utf-16LE", errors="replace")
    with open(
        os.path.join(jobdir, f"output_{job_idx}.log"), "w", encoding="utf-16le"
    ) as f:
        f.write(data)


def execute_olex2c_wrapper(args):
    return execute_olex2c(*args)


if __name__ == "__main__":
    # This script should be started from the folder containing the .cif and ins file!
    base_dir = os.getcwd()

    if len(sys.argv) < 2:
        print(
            "Please provide a folder name, where the calculation should be performed."
        )
        sys.exit(1)

    work_dir = os.path.join(base_dir, sys.argv[1])

    from_the_beginning = True
    # Check if sys.argv[1] is a folder and contains a *.sisy file
    if os.path.isdir(work_dir):
        if any(file.endswith(".sisy") for file in os.listdir(work_dir)):
            print(f"Found a .sisy file in {work_dir}. Continuing run from bevore.")
            from_the_beginning = False
    else:
        print(f"{work_dir} does not exist, creating it.")
        os.makedirs(work_dir)

    n_jobs = 0
    if from_the_beginning:
        n_jobs = write_sisy_file(work_dir, "sisyphos.sisy")
    else:
        # Check if the .sisy file is complete
        with open(os.path.join(work_dir, "sisyphos.sisy"), "r") as f:
            for line in f:
                if line.startswith("#") or line == "\n":
                    continue
                n_jobs += 1

    exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "olexsys-olex2c.exe")

    exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "olexsys-olex2c.exe")

    os.environ["OLEX2_CCTBX_DIR"] = os.path.join(os.path.abspath(__file__), "cctbx")
    os.environ["LOAD_HEADLESS_PLUGINS"] = "True"
    # This Path needs to be adjusted to the location of the ORCA executable
    os.environ["PATH"] = "C:\\ORCA_6.0.0;" + os.environ["PATH"]

    os.environ["SISYPHOS_base_path"] = base_dir
    os.environ["SISYPHOS_work_path"] = work_dir
    # #creating a methodfile list and file

    n_workers = ncpus // cpus_per_job
    os.environ["ncpus"] = str(cpus_per_job)

    print("Number of jobs: ", n_jobs)
    print("Number of workers: ", n_workers)
    print("Number of cpus per job: ", cpus_per_job)
    # Create a list of arguments for each call to execute_olex2c
    args_list = [(exe, work_dir, indexes) for indexes in range(n_jobs)]
    # # Create a multiprocessing Pool with ncpus workers
    with multiprocessing.Pool(n_workers) as pool:
        # Use the pool to map execute_olex2c to the args_list
        pool.map(execute_olex2c_wrapper, args_list)

    print("Done.")
