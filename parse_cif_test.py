import os

elements = ["Au"]

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
              for elem in elements:
                switch = True                
                if line.startswith(f" {elem}  ") & switch:
                    print(line)
                    if " . . " in line:
                      continue
    return out
            
df1 = parse_cif("./pirEt_epsilon_11929_smart3_emp.cif")
print(df1)