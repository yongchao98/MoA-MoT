def solve_reaction():
    """
    Determines the IUPAC name of the product of 1,3-dibromo-2-iodobenzene 
    and excess phenyl magnesium bromide.
    """
    
    # 1. Define the initial state of the substituted benzene ring.
    # The keys represent the positions on the benzene ring.
    # The values represent the substituents.
    substituents = {1: 'Bromo', 2: 'Iodo', 3: 'Bromo'}
    
    print("Starting molecule: 1,3-dibromo-2-iodobenzene")
    print("Reagent: Phenyl magnesium bromide (excess)")
    print("Reaction Conditions: Reflux in THF, followed by aqueous work-up.")
    print("-" * 40)
    print("The reaction proceeds by substituting the halogens with phenyl groups.")
    print("The order of reactivity for the halogens is I > Br.")
    print("-" * 40)

    # 2. First substitution: The most reactive halogen (Iodo) is replaced.
    position_to_replace = 2  # Iodine is at position 2 and is most reactive
    substituents[position_to_replace] = 'Phenyl'
    print("Step 1: The highly reactive Iodo group at position 2 is replaced by a Phenyl group.")
    # The intermediate is 1,3-dibromo-2-phenylbenzene.
    print("Intermediate product: 1,3-dibromo-2-phenylbenzene")
    
    # 3. Second substitution: One of the Bromo groups is replaced.
    position_to_replace = 1  # Bromine at position 1 is replaced next.
    substituents[position_to_replace] = 'Phenyl'
    print("\nStep 2: One Bromo group (at position 1) is replaced by a Phenyl group.")
    # The intermediate is 3-bromo-1,2-diphenylbenzene. IUPAC name is 1-bromo-2,3-diphenylbenzene.
    print("Intermediate product: 1-bromo-2,3-diphenylbenzene")
    
    # 4. Third substitution: The final Bromo group is replaced.
    position_to_replace = 3 # Final bromine at position 3 is replaced.
    substituents[position_to_replace] = 'Phenyl'
    print("\nStep 3: The final Bromo group (at position 3) is replaced by a Phenyl group.")

    # 5. Determine the final product's name
    # The final substituents are Phenyl groups at positions 1, 2, and 3.
    # The name is 1,2,3-triphenylbenzene.
    final_name = "1,2,3-triphenylbenzene"
    
    print("-" * 40)
    print(f"All halogens have been substituted. The final structure has three phenyl groups at positions 1, 2, and 3.")
    print("Final Product IUPAC Name: 1,2,3-triphenylbenzene")

    # This part of the code demonstrates how the final name contains numbers as requested.
    print("\nPrinting each number in the final name:")
    for char in final_name:
        if char.isdigit():
            print(char)

solve_reaction()