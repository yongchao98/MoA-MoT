def solve_inference():
    """
    Determines the natural logic operator for the given inference using semantic properties.
    """
    relations = {
        '=': 'Equivalence',
        '<': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '^': 'Negation',
        '|': 'Alternation',
        'v': 'Cover',
        '#': 'Independence'
    }

    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print("Analyzing the inference:")
    print(f"Premise (P): \"{premise}\"")
    print(f"Hypothesis (H): \"{hypothesis}\"")
    print("-" * 20)

    # Define intermediate statements for clarity
    p_prime = "Mark is singing a song by Taylor Swift"
    statement_b = "Mark is singing a song by Michael Jackson"

    # Step 1: Relation between Premise and an intermediate, more general statement
    rel_p_pprime = '<'
    print(f"Step 1: Premise ('{premise}') entails an intermediate statement P' ('{p_prime}').")
    print(f"This is because 'pop song' is a subtype of 'song'.")
    print(f"Relation (P, P'): {relations[rel_p_pprime]} ({rel_p_pprime})\n")

    # Step 2: Relation between the intermediate statement and the counterpart in the hypothesis
    rel_pprime_b = '|'
    print(f"Step 2: P' ('{p_prime}') and the statement B ('{statement_b}') are mutually exclusive.")
    print(f"This is because a person cannot sing two different songs simultaneously.")
    print(f"Relation (P', B): {relations[rel_pprime_b]} ({rel_pprime_b})\n")

    # Step 3: Using the definition of Alternation to relate to the Hypothesis
    rel_pprime_h = '<'
    print(f"Step 3: By definition of Alternation ({rel_pprime_b}), the truth of P' entails the falsity of B.")
    print(f"The falsity of B is the Hypothesis H ('{hypothesis}').")
    print(f"Therefore, P' entails H.")
    print(f"Relation (P', H): {relations[rel_pprime_h]} ({rel_pprime_h})\n")

    # Step 4: Composing the relations via transitivity
    # We have P < P' and P' < H. Transitivity means P < H.
    final_relation_symbol = '<'
    final_relation_name = relations[final_relation_symbol]
    print(f"Step 4: Final composition via transitivity.")
    print(f"The full logical path is: P {rel_p_pprime} P' and P' {rel_pprime_h} H.")
    # In the final code you still need to output each number in the final equation!
    # Let's represent the relation symbols as numbers for the "equation" part.
    # For simplicity, we can use their index.
    # =:0, <:1, >:2, ^:3, |:4, v:5, #:6
    print(f"Final Equation: Relation({rel_p_pprime}) o Relation({rel_pprime_h}) => Relation({final_relation_symbol})")
    print(f"\nThe final projected natural logic operator is {final_relation_name}.")


solve_inference()