def solve_chemistry_problem():
    """
    Determines the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene with excess phenyl magnesium bromide.
    """

    # 1. Define reactants and conditions
    aryl_halide = {
        "name": "1,3-dibromo-2-iodobenzene",
        "substituents": {
            1: "Bromo",
            2: "Iodo",
            3: "Bromo"
        }
    }
    grignard_reagent = "Phenyl Magnesium Bromide"
    nucleophile = "Phenyl"
    conditions = "Excess reagent and reflux"

    print("--- Reaction Analysis ---")
    print(f"Starting Material: {aryl_halide['name']}")
    print(f"Reagent: {grignard_reagent} ({conditions})")
    print("Reaction Type: Grignard Cross-Coupling")
    print("-" * 25)

    # 2. Analyze halogen reactivity
    print("\n--- Substitution Steps ---")
    print("The reactivity of halogens in this reaction is I > Br > Cl.")
    
    # Simulate the reaction by replacing halogens with the nucleophile
    final_substituents = aryl_halide["substituents"].copy()

    # Step 1: Replace Iodine (most reactive)
    print("Step 1: The 'Iodo' group at position 2 is replaced by a 'Phenyl' group.")
    final_substituents[2] = nucleophile
    
    # Step 2 & 3: Replace Bromines (less reactive, but conditions are harsh)
    print("Step 2: The 'Bromo' group at position 1 is replaced by a 'Phenyl' group.")
    final_substituents[1] = nucleophile
    print("Step 3: The 'Bromo' group at position 3 is replaced by a 'Phenyl' group.")
    final_substituents[3] = nucleophile

    # 3. Determine the final product and its IUPAC name
    parent_molecule = "benzene"
    substituent_name = "phenyl"
    
    # Get sorted positions of the new substituents
    positions = sorted(final_substituents.keys())
    
    # Determine the prefix for the number of substituents (e.g., di, tri, tetra)
    count_prefix_map = {1: "", 2: "di", 3: "tri", 4: "tetra"}
    count = len(positions)
    prefix = count_prefix_map.get(count, "")
    
    # Format the positions as a comma-separated string
    position_str = ",".join(map(str, positions))

    final_iupac_name = f"{position_str}-{prefix}{substituent_name}{parent_molecule}"

    print("\n--- Final Product ---")
    print(f"All halogen atoms have been substituted by '{nucleophile}' groups.")
    print(f"The final product has '{nucleophile}' groups at positions: {position_str}.")
    print("\n-------------------------------------------")
    print(f"The IUPAC name of the final product is: {final_iupac_name}")
    print("The numbers in the name are: 1, 2, and 3.")
    print("-------------------------------------------")


# Run the solver
solve_chemistry_problem()
