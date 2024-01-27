import pandas as pd
from scipy import constants as conts

def ret_energy(wavelength):
    return (conts.h * conts.c) / (wavelength*conts.e) / conts.angstrom

def ret_wl(energy):
    return (conts.h * conts.c) / (energy*conts.e) / conts.angstrom

def evaluate(loc):
    with open(loc, "r") as inp:
        df = pd.DataFrame()
        blocks = []
        buffer = []
        for line in inp:
            if line.startswith("++++"):
                blocks.append(buffer)
                buffer = []
            else:
                buffer.append(line)
    for block in blocks:
        outdir = {}
        for i in range(len(block)):
            if "\t" in block[i]:
                block[i] = block[i].split("\t")[1]
                block[i] = block[i].strip("\n")
        for elem in block:
            for piece in elem.split(";"):
                if ":" in piece:
                    key, value = piece.split(":")
                    outdir[key] = value
                else:
                    continue
        df = df.append(outdir, ignore_index =True)
    df["energy"] = ret_energy(df["wavelength"])
    df.to_csv("SYSout.csv")
