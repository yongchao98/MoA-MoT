def solve_natural_logic():
    """
    Solves for the final natural logic operator based on MacCartney's framework.
    """
    # Define relation symbols for clarity in the output
    relations = {
        "Forward Entailment": "⊂",
        "Alternation": "∧"
    }

    # Define the premise and hypothesis
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print("Step-by-step derivation of the natural logic relation:")
    print("-" * 50)
    print(f"Premise (P): \"{premise}\"")
    print(f"Hypothesis (H): \"{hypothesis}\"")
    print("-" * 50)

    # Step 1: Analyze the hypothesis without the negation
    print("\nStep 1: Determine the relation between the Premise (P) and the positive form of the Hypothesis (H').")
    hypothesis_positive = "Mark is singing a song by Michael Jackson"
    print(f"Let H' = \"{hypothesis_positive}\"\n")

    # Decompose the sentences to find the edits
    edit1_p = "a pop song"
    edit1_h = "a song"
    rel1_name = "Forward Entailment"
    rel1_sym = relations[rel1_name]

    edit2_p = "Taylor Swift"
    edit2_h = "Michael Jackson"
    rel2_name = "Alternation"
    rel2_sym = relations[rel2_name]

    print("To transform P into H', two main edits are needed:")
    print(f"  1. Edit '{edit1_p}' → '{edit1_h}'. The relation is {rel1_name} ({rel1_sym}).")
    print(f"  2. Edit '{edit2_p}' → '{edit2_h}'. The relation is {rel2_name} ({rel2_sym}).\n")

    # Composition of the edits
    intermediate_rel_name = "Alternation"
    intermediate_rel_sym = relations[intermediate_rel_name]
    print("Composing these two edits using MacCartney's join table:")
    print(f"  Final Equation for Step 1: join({rel1_sym}, {rel2_sym}) = {intermediate_rel_sym}")
    print(f"Result: The relation between P and H' is {intermediate_rel_name} ({intermediate_rel_sym}).\n")
    print("-" * 50)

    # Step 2: Account for the negation in the original hypothesis
    print("\nStep 2: Account for the negation in the original Hypothesis H.")
    print(f"The final Hypothesis H is the negation of H', i.e., H = not H'.")
    print(f"We need to find the relation between P and (not H').\n")

    print(f"An {intermediate_rel_name} ({intermediate_rel_sym}) relation means P and H' are mutually exclusive.")
    print("Therefore, if P is true, H' must be false. This means P entails (not H').")
    print("This one-way implication defines the Forward Entailment relation.\n")

    # Final relation
    final_rel_name = "Forward Entailment"
    final_rel_sym = relations[final_rel_name]
    print("The final projected natural logic operator is therefore " + final_rel_name + ".")

    print("\nThe final equation is:")
    print(f"  relation(P, H) = neg(relation(P, H')) = neg({intermediate_rel_sym}) = {final_rel_sym}")

solve_natural_logic()