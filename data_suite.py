"""
Currently not used!!!
Project to be loaded as a module into custom scripts instead.
"""

raise NotImplementedError("This script is currently not used! Load this package into your own script with 'import SIF'")

###########################################
# ~~~ Simulation Input Finagler ~~~       #
# Tools to mess around                    #
# with assorted LAMMPS data               #
#                                         #
# Written by Kieran Harris                #
###########################################

# Import each tool that can be run (from ./tools/)
from tools import *

# List of name & function for each tool
tools = [
    {"name": "Data->template conversion",
     "function": data_to_template.main},

    {"name": "Template charge fixing",
     "function": None}
    ]


def main():
    
    print("Welcome to the (unofficial) LAMMPS data suite!")

    # Keep asking until we get a valid answer and break
    while True:
        print("")
        print("Choose which tool you would like to use:")
        # Iterate over each tool for options menu
        for i,tool in enumerate(tools):
            print(f"{i} - {tool['name']}")

        # Get option input & check it's valid
        try:
            opt = int(input(""))
            if 0 <= opt < len(tools):
                break
            else:
                print("\nPlease input a valid integer.")
        except ValueError:  # Integer conversion failed
            print("\nPlease input an integer.")

    print(f"\nChose {tools[opt]['name']}")

    # Run main function of chosen tool
    tools[opt]["function"]()

    exit()

main()
