def solve_reaction():
    """
    This function determines and prints the product of a two-step organic reaction.
    """
    # Define the reaction components
    starting_material = "N,N-diethyl-3-dimethylaminobenzamide"
    reagents_step_1 = "1) sec-BuLi, TMEDA, THF"
    reagents_step_2 = "2) Methyl iodide (CH3I)"

    # Explain the reaction mechanism
    print("Reaction Analysis:")
    print("The reaction is a Directed ortho-Metalation (DoM) followed by an electrophilic quench.")
    print("Step 1: The strong base, sec-BuLi, selectively removes a proton from the benzene ring.")
    print("The deprotonation occurs at position 2, which is 'ortho' to both the N,N-diethylamide and the 3-dimethylamino directing groups.")
    print("Step 2: The resulting nucleophilic intermediate attacks the electrophile, methyl iodide.")
    
    # Define the final product based on the mechanism
    # The result is the addition of a methyl group at position 2.
    final_product = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"

    # Print the final result
    print("\n---")
    print("Starting Material: {}".format(starting_material))
    print("Final Product: {}".format(final_product))

solve_reaction()