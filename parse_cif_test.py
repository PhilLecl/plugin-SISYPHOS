import os

elem_string = ["U"]

dat_names = ["mu", 
             "wavelength", 
             "F000", 
             "tot_reflIns", 
             "goof", 
             "R_all", 
             "R1", 
             "wR2", 
             "last Shift"]

corr_filts = ["exptl_absorpt_coefficient_mu", 
              "diffrn_radiation_wavelength", 
              "exptl_crystal_F_000",
              "diffrn_reflns_number",
              "refine_ls_goodness_of_fit_ref",
              "refine_ls_R_factor_all",
              "refine_ls_R_factor_gt",
              "refine_ls_wR_factor_ref",
              "REM Shift_max"]

def parse_cif(loc):
    out = {}
    with open(loc, "r") as incif:
        for line in incif:
            for i,filter in enumerate(corr_filts):
                if filter in line:
                    out[f"{dat_names[i]}"] = float(line.split()[-1])

    with open(loc, "r") as incif:
        switch2 = False
        for line in incif:
            for elem in elem_string:
                switch = True                
                if  line.startswith(f" {elem} ") & switch:
                    switch = False
                    if "(" in line.split(" ")[3]:
                        print(line)
                        fp = (float(line.split(" ")[2].split("(")[0]), int(line.split(" ")[2].split("(")[1][:-1]))
                        fdp = (float(line.split(" ")[3].split("(")[0]), int(line.split(" ")[3].split("(")[1][:-1]))
                    else:
                        print(line)
                        fp = float(line.split(" ")[2])
                        fdp = float(line.split(" ")[3])
                    out[f"{elem}_anoms"] = (fp,fdp)

            if line.startswith("  _atom_site_refinement_flags_occupancy"):
                switch2 = True
                continue
            if switch2:
                if line.startswith("\n"):
                    switch2 = False
                else:
                    lin = line.split(" ")
                    atom = lin[1]
                    ueq = lin[6].split("(")[0]
                    ueq_delta = lin[6].split("(")[1][:-1]
                    out[f"{atom}_ueq"] = (float(ueq), int(ueq_delta))
    return out
            
df1 = parse_cif("./OW3_rt_Au5filt_17280_run.cif")
print(df1)