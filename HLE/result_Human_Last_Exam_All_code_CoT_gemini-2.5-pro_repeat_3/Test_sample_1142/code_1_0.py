def get_electron_configuration(atomic_number):
    """
    Determines the electron configuration based on the Aufbau principle.
    Returns the configuration as a dictionary and identifies the
    outermost, partially-filled subshell.
    """
    # Standard orbital filling order and their maximum electron capacities
    orbitals = [
        ('1s', 2), ('2s', 2), ('2p', 6), ('3s', 2), ('3p', 6),
        ('4s', 2), ('3d', 10), ('4p', 6), ('5s', 2), ('4d', 10), ('5p', 6)
    ]
    
    config = {}
    electrons_left = atomic_number
    last_partially_filled_orbital = None
    
    for name, capacity in orbitals:
        if electrons_left == 0:
            break
        
        electrons_to_fill = min(electrons_left, capacity)
        config[name] = electrons_to_fill
        electrons_left -= electrons_to_fill
        
        # This is the highest energy subshell being filled
        if electrons_to_fill > 0:
            last_partially_filled_orbital = (name, capacity, electrons_to_fill)
            
    return config, last_partially_filled_orbital

def calculate_valence(orbital_info):
    """
    Calculates valence based on the rule of completing the subshell.
    """
    name, capacity, current_electrons = orbital_info
    # Valence is the number of electrons needed to complete the subshell
    valence = capacity - current_electrons
    return valence

def main():
    """
    Main function to execute the analysis for NiC.
    """
    # Define atomic numbers
    atomic_numbers = {'C': 6, 'Ni': 28}

    valences = {}

    # --- Step 1 & 2: Calculate Configurations and Valences ---
    print("--- Step 1 & 2: Calculating Valences ---")
    for symbol, z in atomic_numbers.items():
        config, orbital_info = get_electron_configuration(z)
        valence = calculate_valence(orbital_info)
        valences[symbol] = valence
        
        name, capacity, current = orbital_info
        print(f"Analysis for {symbol} (Z={z}):")
        print(f"  Electron Configuration: {config}")
        print(f"  Outermost partially filled subshell: {name}")
        print(f"  This subshell has {current} of {capacity} electrons.")
        print(f"  Valence = electrons to complete subshell = {capacity} - {current} = {valence}")
        print("-" * 20)

    # --- Step 3: Analyze NiC Bonding ---
    print("--- Step 3: Analyzing Bonding in NiC ---")
    valence_C = valences['C']
    valence_Ni = valences['Ni']
    average_valence = (valence_C + valence_Ni) / 2
    
    print(f"Covalent bonding requirements: C wants {valence_C} bonds, Ni wants {valence_Ni} bonds.")
    print("For a 1:1 NiC crystal, the atoms likely form a compromise structure.")
    print("The average number of bonds (coordination number) per atom is:")
    print(f"  ({valence_C} + {valence_Ni}) / 2 = {average_valence:.1f}")
    print("-" * 20)

    # --- Step 4: Select Crystal Structure ---
    print("--- Step 4: Selecting the Crystal Structure ---")
    print(f"The average coordination number of {int(average_valence)} points to a structure where all atoms form {int(average_valence)} bonds.")
    print("Reviewing the options:")
    print("  A. flattened tetrahedral (4)")
    print("  B. tiling by squares (4)")
    print("  C. tiling by octagons and squares (3)")
    print("  D. tiling by hexagons (3)")
    print("  E. foliation by chains (2)")
    print("  F. partition into rings (2)")
    print("\nOptions C and D have the required coordination number of 3.")
    print("A 3-coordinate structure is most stable with trigonal planar geometry (120-degree bond angles).")
    print("Option D (hexagonal tiling) provides perfect 120-degree angles.")
    print("Option C (octagons and squares) has strained angles of 90 and 135 degrees.")
    print("Therefore, the hexagonal tiling (D) is the most likely structure.")
    print("-" * 20)

    # --- Step 5: Assess Shear Strength ---
    print("--- Step 5: Assessing Shear Strength Isotropy ---")
    print("The predicted structure is a hexagonal lattice.")
    print("A hexagonal lattice is fundamentally anisotropic; its mechanical properties, including shear strength, are different along different directions (e.g., 'armchair' vs. 'zigzag').")
    print("Therefore, the crystal shear strength is not isotropic.")
    print("The answer to the second question is 'no'.")
    print("-" * 20)

    print("Final Answer Summary:")
    print("Choice: D (tiling by hexagons)")
    print("Isotropic?: no")

main()