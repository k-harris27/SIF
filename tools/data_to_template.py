import shutil
from utils import *


def main():
    
    print("\nWARNING: DEPRECATED")
    print("A new version should be in the works...")
    
    print("\nEnter species name (input file should be species_frag.data):")
    species=input("")
    
    # Create temp file we can work on
    shutil.copyfile(species+'_frag.data','temp.frag')
    
    # Open temp file
    with open('temp.frag',) as data_file:
        # Read the file line by line into a list to edit
        data = data_file.readlines()

    remove_words = ['types','xlo xhi','ylo yhi','zlo zhi']
    remove_sections = ['Masses','Velocities']
    
    # Remove all lines containing words in remove_words
    data = [line for line in data if not any(match in line for match in remove_words)]
    
    # Now iterate over each line in the list.
    
    # template file output lines in list format
    t_data=["Template file for reactive part of species {}\n".format(species),"\n"]
    # Atom data will need post-processing, so we set it to the side
    atom_data=[]
    remove=False  # Are we in a section to be removed?
    atoms=False  # Are we in the atoms section (which needs special work)?
    # Remove sections we don't want. Skip line 0, which is just comments.
    for line in data[1:]:
        # Remove any comments (after a #)
        line,_,_ = line.partition('#')
        
        line = line.strip()  # Remove newline and leading/trailing whitespace
    
        if not line: continue  # If line is blank, skip it.
        
        # Section headers are the only lines to start with a capital letter.
        # If we are in a new section, reset removal flags etc
        if line[0].isupper():
    
            # If we have reached a new section and atoms is true
            # then we have reached the end of atoms section.
            if atoms:
                ProcessAtomData(t_data,atom_data)
            remove=False
            atoms=False
    
            if line in remove_sections:
                remove=True
                continue
            elif line == "Atoms":
                atoms=True
                continue
    
            # Newlines go before & after for nicer formatting
            t_data.extend(['\n', line + '\n', '\n'])
            continue
            
        if remove:
            continue
        elif atoms:
            atom_data.append(line)  # Set aside for post-processing
        else:
            t_data.append(line + '\n') # Not header or modified -> keep
    
    # Finally, we write the new template data to file!
    with open(species+"_frag.template","w") as out:
        out.writelines(t_data)
        print("Done!")

    

# Run this upon reaching the end of the Atoms section
def ProcessAtomData(t_data,atom_data):
    coord_lines = []
    type_lines = []
    charge_lines = []
    for atom in atom_data:
        # Atoms are stored in .data files as:
        # ID molID Type charge x y z ix iy iz
        coords = ["","",""]
        images = ["","",""]  # NOTE: Images ignored rn. Coords shouldn't matter.
        # Split by whitespace, assign to relevant variables
        at_id,_,at_type,at_charge,coords[0],coords[1],coords[2],images[0],images[1],images[2] = atom.split()
        coord_lines.append(f"{at_id} {coords[0]} {coords[1]} {coords[2]}\n")
        type_lines.append(f"{at_id} {at_type}\n")
        charge_lines.append(f"{at_id} {at_charge}\n")

    # Append coord section to t_data.
    t_data.extend(["\n","Coords\n","\n"]+coord_lines)
    
    # Append type section to t_data.
    t_data.extend(["\n","Types\n","\n"]+type_lines)

    # Append charge section to t_data.
    t_data.extend(["\n","Charges\n","\n"]+charge_lines)

