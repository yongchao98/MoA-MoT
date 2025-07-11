def solve_reaction_and_name_product():
    """
    This script determines the product of the reaction between
    1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide,
    and then generates its IUPAC name.
    """

    # 1. Define the initial molecule as substituents on a benzene ring
    # The keys are positions, and values are the substituent groups.
    molecule = {1: 'Bromo', 2: 'Iodo', 3: 'Bromo', 4: 'H', 5: 'H', 6: 'H'}

    # 2. Define the relative reactivity of leaving groups (lower number = more reactive)
    # In Grignard cross-coupling, reactivity is I > Br > Cl.
    leaving_group_reactivity = {'Iodo': 1, 'Bromo': 2, 'Chloro': 3}

    # The incoming group from phenyl magnesium bromide
    incoming_group = 'Phenyl'

    print("--- Reaction Analysis ---")
    print(f"Starting Material: 1,3-dibromo-2-iodobenzene")
    print(f"Reagent: Excess Phenyl Magnesium Bromide (source of '{incoming_group}' groups)")
    print("-" * 25)

    # 3. Simulate the reaction step-by-step
    # Create a list of halogens to be replaced, with their positions and reactivity
    halogens_to_replace = []
    for position, group in molecule.items():
        if group in leaving_group_reactivity:
            reactivity = leaving_group_reactivity[group]
            halogens_to_replace.append({'position': position, 'group': group, 'reactivity': reactivity})

    # Sort the list by reactivity to ensure the most reactive group is replaced first
    halogens_to_replace.sort(key=lambda x: x['reactivity'])

    # Perform the substitutions
    step = 1
    for halogen in halogens_to_replace:
        pos = halogen['position']
        group = halogen['group']
        print(f"Step {step}: The most reactive leaving group is '{group}' at position {pos}.")
        print(f"           '{incoming_group}' group substitutes '{group}'.")
        molecule[pos] = incoming_group
        step += 1

    print("-" * 25)
    print("All halogens have been substituted due to excess reagent.")
    
    # 4. Generate the IUPAC name for the final product
    
    substituent_positions = []
    substituent_name = ""
    for position, group in molecule.items():
        if group == incoming_group:
            substituent_positions.append(str(position))
            substituent_name = group.lower()

    count_prefix_map = {1: "", 2: "di", 3: "tri", 4: "tetra"}
    count = len(substituent_positions)
    prefix = count_prefix_map.get(count, "poly")

    # Final name construction
    locants = ",".join(substituent_positions)
    # Final IUPAC name = [locants]-[prefix][substituent_name][parent_chain]
    final_name = f"{locants}-{prefix}{substituent_name}benzene"
    
    print(f"Final Product Structure: Benzene ring with '{incoming_group}' groups at positions {locants}.")
    print("\n--- Final IUPAC Name ---")
    print("The final name is constructed as: [positions]-[prefix][substituent]benzene")
    print(f"Positions: {substituent_positions[0]}, {substituent_positions[1]}, {substituent_positions[2]}")
    print(f"Number of groups: {count} -> Prefix: {prefix}")
    print(f"Substituent: {substituent_name}")
    print(f"Parent: benzene")
    print("\nResulting Name:")
    print(final_name)


solve_reaction_and_name_product()