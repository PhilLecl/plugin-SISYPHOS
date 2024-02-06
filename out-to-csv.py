import pandas as pd
import sys
from scipy import constants as conts

def ret_energy(wavelength):
    return (conts.h * conts.c) / (wavelength*conts.e) / conts.angstrom

def ret_wl(energy):
    return (conts.h * conts.c) / (energy*conts.e) / conts.angstrom

def main():
    for arg in sys.argv:
        evaluate(arg)

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
            if "\n" in block[i]:
                block[i] = block[i].strip("\n")
        for elem in block:
            for piece in elem.split(";"):
                if ":" in piece:
                    key, value = piece.split(":")
                    outdir[key] = value
                else:
                    continue
        result = pd.DataFrame(outdir, index=[0])
        df = pd.concat([df, result], ignore_index=True)
    df.to_csv("SYSout.csv")

if __name__ == '__main__':
    main()