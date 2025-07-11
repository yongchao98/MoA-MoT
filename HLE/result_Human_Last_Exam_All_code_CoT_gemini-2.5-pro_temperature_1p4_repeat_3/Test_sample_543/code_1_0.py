def solve_grignard_reaction():
    """
    Determines the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene with excess phenyl magnesium bromide.
    """

    # Step 1: Define the starting molecule's substituents and the incoming group.
    # We represent the substituted benzene ring as a dictionary for clarity.
    # Positions are keys, and substituent names are values.
    substituents = {1: 'Bromo', 2: 'Iodo', 3: 'Bromo'}
    reagent_group = 'Phenyl'

    print("Starting material: 1,3-dibromo-2-iodobenzene")
    print(f"Reagent: Excess Phenyl magnesium bromide (adds '{reagent_group}' groups)")
    print("Conditions: Reflux in THF, followed by aqueous work-up.")
    print("-" * 50)

    # Step 2: Define the reactivity order. The C-I bond is the most reactive,
    # followed by the C-Br bonds.
    # The script will replace halogens in this specific order.
    
    print("Reaction proceeds based on halogen reactivity (I > Br).")
    print("Since the reagent is in excess, all halogens will be substituted.\n")

    # Step 3: Simulate the step-by-step substitution reaction.
    
    # First, replace the most reactive halogen: Iodine.
    iodo_position = [pos for pos, sub in substituents.items() if sub == 'Iodo'][0]
    print(f"Step 1: The 'Iodo' group at position {iodo_position} is the most reactive and is replaced first.")
    substituents[iodo_position] = reagent_group
    
    # Second, replace the remaining halogens: Bromine.
    bromo_positions = sorted([pos for pos, sub in substituents.items() if sub == 'Bromo'])
    print(f"Step 2: The 'Bromo' group at position {bromo_positions[0]} is replaced.")
    substituents[bromo_positions[0]] = reagent_group
    
    print(f"Step 3: The final 'Bromo' group at position {bromo_positions[1]} is replaced.")
    substituents[bromo_positions[1]] = reagent_group
    
    print("-" * 50)

    # Step 4: Construct the IUPAC name for the final product.
    print("All halogen substituents have been replaced by phenyl groups.")
    
    # Get the sorted positions of the final phenyl groups.
    phenyl_positions = sorted(substituents.keys())

    # Determine the multiplier prefix for the number of substituents.
    count = len(phenyl_positions)
    prefixes = {1: '', 2: 'di', 3: 'tri', 4: 'tetra'}
    prefix = prefixes.get(count, 'poly')

    # Format the position numbers (locants).
    locants_str = ",".join(map(str, phenyl_positions))

    # Assemble the final name.
    base_name = "phenylbenzene"
    final_name = f"{locants_str}-{prefix}{base_name}"
    
    print("\n--- Final Product ---")
    print("The final molecule is a benzene ring with three phenyl groups.")
    print(f"The numbers (locants) indicating the positions of the phenyl groups are: {phenyl_positions[0]}, {phenyl_positions[1]}, and {phenyl_positions[2]}.")
    print("\nFinal IUPAC Name:")
    print(final_name)


if __name__ == "__main__":
    solve_grignard_reaction()