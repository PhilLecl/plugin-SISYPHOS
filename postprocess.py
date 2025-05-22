import glob
import os

import pandas as pd

# results.txt looks like this (without the #):
# NoSpherA2_Dict:
# basis_name:def2-TZVP,method:r2SCAN,ncpus:4,mem:2,charge:0,multiplicity:0,full_HAR:True,Max_HAR_Cycles:15,becke_accuracy:Normal,Relativistic:False,h_aniso:True,h_afix:False,add_disp:None,cluster_radius:0,DIIS:0.01,cluster_grow:True,ORCA_SCF_Conv:NoSpherA2SCF,ORCA_SCF_Strategy:EasyConv,ORCA_Solvation:Water,pySCF_Damping:0.6,ORCA_DAMP:False,NoSpherA2_RI_FIT:True,auxiliary_basis:auto_aux,auto_aux_beta:1.7,
# Stats-GetHklStat:
# TotalReflections:17517,UniqueReflections:3190,DataCount:3190,FriedelOppositesMerged:0,InconsistentEquivalents:4,SystematicAbsencesRemoved:0,MinD:0.5838365560882717,MaxD:8.616100000000001,LimDmin:0.355365,LimDmax:100.0,FilteredOff:0,OmittedByUser:0,OmittedReflections:0,IntensityTransformed:0,Rint:0.02853929076435011,Rsigma:0.030547274979897555,MeanIOverSigma:28.156968919312433,Completeness:0.9962523422860712,MaxIndices:(8, 11, 29),MinIndices:(-8, 0, 0),FileMaxIndices:(8, 10, 29),FileMinIndices:(-8, -11, -29),ReflectionAPotMax:24,FriedelPairCount:2076,Redundancy:(3485, 2250, 1276, 631, 304, 125, 59, 19, 7, 13, 7, 4, 1, 1),
# Cell-Stats:
# a:5.1288(3),b:6.8814(3),c:17.2322(11),alpha:90,beta:90,gamma:90,volume:608.18(6),Z:4,Zprime:1,
# CIF-stats:
# mu:0.14,F000:312.28,tot_reflIns:17517.0,wavelength:0.71073,goof:1.0681,R_all:0.0365,R1:0.0271,wR2:0.0466,last Shift:-0.0006,
# refine_dict:
# max_peak:0.2246993116959985,max_hole:-0.203670959276601,res_rms:0.05742005608646172,goof:1.0681,max_shift_over_esd:-0.0005670637862517267,hooft_str:-0.2(3),R1_all:0.0365313469967893,R1_gt:0.027074484617873998,wR2:0.046627292199113456,cycles:15,time:183.0656977,
# bondlengths:
# O001-C006:1.2603907805960624,O002-H002:1.010433718643144,O002-C007:1.3111519250732477,O003-C006:1.2469978762953324,O004-C007:1.2221357609806656,N005-H00c:0.9974267560154373,N005-H00a:1.0262564704112738,N005-H00b:1.026936104231928,N005-C009:1.4888211183897881,C006-C009:1.5328733868826363,C007-C008:1.5052740773491122,C008-H00d:1.0806746155314837,C008-H00e:1.091423746151058,C008-C00A:1.524178261417927,C009-H009:1.0692821303189013,C009-C00A:1.5300361012785306,C00A-H00f:1.0733157337164148,C00A-H00g:1.0750455638083205,
# bonderrors:
# O001-C006:0.0007379390584126034,O002-H002:0.013579715575292074,O002-C007:0.0007731842347214881,O003-C006:0.0007577036964648571,O004-C007:0.0007857163677244854,N005-H00c:0.012230506570043907,N005-H00a:0.012021914131040392,N005-H00b:0.010189082059966886,N005-C009:0.0008133693851701788,C006-C009:0.0008256568592481016,C007-C008:0.0008755313270173721,C008-H00d:0.010118728777790423,C008-H00e:0.010001853575319648,C008-C00A:0.0009217271196228342,C009-H009:0.009293642484198618,C009-C00A:0.0008498029270465262,C00A-H00f:0.01039405211845793,C00A-H00g:0.009076795679949835,
# Weight:True
# Nr. NPD:1


class PostProcess:
    def __init__(self, work_dir):
        self.work_dir = work_dir

    def read_refinement_results(self, path_to_result) -> dict:
        """
        Read the refinement results from the specified directory.
        """
        # Check if the file exists
        if not os.path.exists(path_to_result):
            raise FileNotFoundError(
                f"Refinement results file not found: {path_to_result}"
            )

        # Read the file
        with open(path_to_result, "r") as file:
            lines = file.readlines()

        results = {}
        for line in lines:
            if line.endswith(":\n") or "+++++" in line:
                # This is a header line, skip it
                continue
            # Split the line into key-value pairs
            pairs = line.strip().split(",")
            for pair in pairs:
                if ":" in pair:
                    key, value = pair.split(":", 1)
                    key = key.strip()
                    value = value.strip()
                    # Convert numeric values to float if possible
                    try:
                        value = float(value)
                    except ValueError:
                        value = value.strip()
                    results[key] = value
        return results

    def collect_results_to_dataframe(self):
        """Write the results to a pandas DataFrame."""
        # Initialize an empty list to store the results
        all_results = []

        for folder in sorted(
            glob.glob(os.path.join(self.work_dir, "job_*")),
            key=lambda x: int(x.split("_")[-1]),
        ):
            # Read the results from each file
            results = self.read_refinement_results(os.path.join(folder, "results.txt"))
            # Append the results to the list
            all_results.append(results)

        # print(all_results)
        # Create a DataFrame from the list of results
        df = pd.DataFrame(all_results)
        # Save the DataFrame to a CSV file
        output_csv = os.path.join(self.work_dir, "refinement_results.csv")
        df.to_csv(output_csv, index=False)

        print(f"Results saved to {output_csv}")


pp = PostProcess(r"C:\Users\Lukas\Desktop\DGK_BILDER\TEST2\SISY_TEST")
pp.collect_results_to_dataframe()
