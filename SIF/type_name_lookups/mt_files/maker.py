def main():

    with open("../oplsaa_type_to_topo.json", "w") as file_out:
        with open("type_to_topo.lt", "r") as t2t_in:
            file_out.write('{ "type_to_topo" : {\n')
            newline=""
            for line in t2t_in:
                words = line.split()
                type_id = int(words[1][6:])
                topo_ids = words[2].split("_")[1:]
                topo_ids = [id[1:] for id in topo_ids]  # Cut off topo kind identifier
                if len(set(topo_ids)) != 1: raise ValueError("One type with multiple topo ids isn't currently handled!")
                topo_id = topo_ids[0]  # They should all be identical anyway.
                file_out.write(f'{newline}\t"t{type_id:03d}"\t: "tt{topo_id}"')
                newline = ",\n"
            file_out.write("\n},\n")

        with open("dihedrals.lt") as dih_in:
            file_out.write('"dihedrals" : [\n')
            newline=""
            for line in dih_in:
                line = line.split()[0]
                topo_name = line.split(":")[1]
                t_n_str = [f'"tt{n}"' for n in topo_name.split("_")]
                file_out.write(f'{newline}\t[{", ".join(t_n_str)}]')
                newline = ",\n"
            file_out.write("\n]}")

main()