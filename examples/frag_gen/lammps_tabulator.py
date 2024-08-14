import os
file_path = f"{os.path.dirname(__file__)}/ddm_frag.data"

atom_names = {
    1 : "CT",
    2 : "HC",
    3 : "CA",
    4 : "HA",
    5 : "OH",
    6 : "HO",
    7 : "OS",
    8 : "NH2",
    9 : "NHR",
    10 : "NR2",
    11 : "H2N",
    12 : "HNR",
}

def main():
    print("starting")
    with open(file_path) as infile:
        in_atoms = False
        for line in infile:
            line = line.strip()
            if not in_atoms:
                if line == "Atoms":
                    in_atoms = True
                    print("in atoms")
            else:
                words = line.split()
                if len(words) > 1:  # Only true on atom info lines
                    print(f"{int(words[0])-1} & {atom_names[int(words[2])]} & {words[3]} \\\\")
                elif len(words) == 1:
                    break  # Reached the end of the atoms section - done




if __name__ == "__main__":
    main()