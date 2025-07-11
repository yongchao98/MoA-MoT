def solve_natural_logic_inference():
    """
    Solves the entailment problem using MacCartney's compositional framework.
    """

    # The 7 natural logic relations with their symbols and names
    relations = {
        '≡': "Equivalence",
        '<': "Reverse Entailment",
        '>': "Forward Entailment",
        '^': "Negation",
        '|': "Alternation",
        '#': "Cover",
        '~': "Independence"
    }

    # MacCartney's composition table: comp_table[row_relation][column_relation]
    # The order of relations for rows/cols is: ≡, <, >, ^, |, #, ~
    comp_table = {
        '≡': {'≡': '≡', '<': '<', '>': '>', '^': '^', '|': '|', '#': '#', '~': '~'},
        '<': {'≡': '<', '<': '<', '>': '~', '^': '|', '|': '|', '#': '~', '~': '~'},
        '>': {'≡': '>', '<': '~', '>': '>', '^': '#', '|': '~', '#': '#', '~': '~'},
        '^': {'≡': '^', '<': '#', '>': '|', '^': '≡', '|': '>', '#': '<', '~': '~'},
        '|': {'≡': '|', '<': '~', '>': '|', '^': '<', '|': '|', '#': '~', '~': '~'},
        '#': {'≡': '#', '<': '#', '>': '~', '^': '>', '|': '~', '#': '>', '~': '~'},
        '~': {'≡': '~', '<': '~', '>': '~', '^': '~', '|': '~', '#': '~', '~': '~'}
    }

    # Premise and Hypothesis
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print("This script determines the natural logic relation for the given inference.")
    print(f"Premise (P):      \"{premise}\"")
    print(f"Hypothesis (H):   \"{hypothesis}\"")
    print("-" * 20)
    
    print("Plan: We will use a compositional proof by introducing an intermediate sentence (H1).")
    print("The final relation will be the composition of rel(P, H1) and rel(H1, H).")
    print("-" * 20)

    # Define the intermediate sentence H1
    h1 = "Mark is singing a song by Michael Jackson"

    # Step 1: Find the relation between the Premise (P) and H1
    print("Step 1: Find rel(P, H1)")
    print(f"P:  \"{premise}\"")
    print(f"H1: \"{h1}\"")
    print("The semantic relationship is defined by 'a pop song by Taylor Swift' and 'a song by Michael Jackson'.")
    print("These two phrases describe mutually exclusive sets. They are also not exhaustive.")
    rel_p_h1_symbol = '#'
    rel_p_h1_name = relations[rel_p_h1_symbol]
    print(f"The relation for mutually exclusive but non-exhaustive items is '{rel_p_h1_name}'.")
    print(f"rel(P, H1) = {rel_p_h1_name} ({rel_p_h1_symbol})")
    print("-" * 20)

    # Step 2: Find the relation between H1 and the Hypothesis (H)
    print("Step 2: Find rel(H1, H)")
    print(f"H1: \"{h1}\"")
    print(f"H:  \"{hypothesis}\"")
    print("The hypothesis H is the direct logical negation of the sentence H1.")
    rel_h1_h_symbol = '^'
    rel_h1_h_name = relations[rel_h1_h_symbol]
    print(f"The relation between a sentence and its negation is '{rel_h1_h_name}'.")
    print(f"rel(H1, H) = {rel_h1_h_name} ({rel_h1_h_symbol})")
    print("-" * 20)

    # Step 3: Compose the relations to find the final result
    print("Step 3: Compose the relations to find rel(P, H).")
    final_rel_symbol = comp_table[rel_p_h1_symbol][rel_h1_h_symbol]
    final_rel_name = relations[final_rel_symbol]
    print(f"The final equation is comp(rel(P, H1), rel(H1, H))")
    print(f"Substituting the symbols: comp({rel_p_h1_symbol}, {rel_h1_h_symbol}) = {final_rel_symbol}")
    print(f"Substituting the names:   comp({rel_p_h1_name}, {rel_h1_h_name}) = {final_rel_name}")
    print("-" * 20)
    
    print(f"The name of the final projected natural logic operator is: {final_rel_name}")

solve_natural_logic_inference()