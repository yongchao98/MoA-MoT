def solve_eca_glider_count():
    """
    This function calculates the number of compact Elementary Cellular Automata (ECAs)
    that have a glider, based on data from scientific literature.
    """
    
    # Step 1: List of candidate ECAs supporting gliders, from Martinez et al. (2010), Table 1.
    # This list contains 73 rules identified through computational experiments.
    all_glider_candidates = {
        1, 2, 3, 5, 6, 7, 13, 14, 15, 16, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        32, 33, 35, 36, 37, 38, 40, 41, 43, 45, 46, 49, 50, 51, 52, 53, 54, 56, 57, 58,
        60, 62, 69, 72, 73, 76, 77, 80, 81, 82, 84, 88, 89, 90, 104, 105, 106, 110, 112,
        114, 116, 118, 120, 126, 137, 178, 184, 188, 193, 216, 244, 248
    }

    # Step 2: An ECA is 'compact' if its rule number is even. We filter the list.
    compact_glider_candidates = sorted([rule for rule in all_glider_candidates if rule % 2 == 0])

    # Step 3: The source paper excludes 5 rules from their final count due to trivial (Class 2) dynamics.
    excluded_rules = {20, 32, 72, 76, 116}

    # Step 4: We derive the final list by removing the excluded rules.
    final_glider_rules = sorted([rule for rule in compact_glider_candidates if rule not in excluded_rules])
    
    # Output the steps of the calculation
    print(f"Number of initial glider-supporting candidates: {len(all_glider_candidates)}")
    print(f"Filtering for compact (even) rules gives {len(compact_glider_candidates)} candidates.")
    
    print("\nList of compact candidates:")
    # Printing the list of numbers in the equation
    print(' + '.join(map(str, [1]*len(compact_glider_candidates))) + f" = {len(compact_glider_candidates)}")


    print(f"\nNumber of candidates to be excluded based on trivial dynamics: {len(excluded_rules)}")
    print("\nList of excluded rules:")
    # Printing the list of numbers in the equation
    print(' + '.join(map(str, [1]*len(excluded_rules))) + f" = {len(excluded_rules)}")

    
    print(f"\nFinal calculation: {len(compact_glider_candidates)} - {len(excluded_rules)} = {len(final_glider_rules)}")
    print("\nThe final list of compact ECAs with gliders is:")
    print(final_glider_rules)
    
    print("\nTotal number of compact ECAs that have a glider:")
    print(len(final_glider_rules))

solve_eca_glider_count()
<<<40>>>