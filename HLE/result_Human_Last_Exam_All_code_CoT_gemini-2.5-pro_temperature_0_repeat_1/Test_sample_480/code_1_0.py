def solve_entailment():
    """
    This function explains the step-by-step derivation of the natural logic operator
    for the given premise and hypothesis, following MacCartney's framework.
    """
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"
    unnegated_hypothesis = "Mark is singing a song by Michael Jackson"

    print(f"Premise (P): \"{premise}\"")
    print(f"Hypothesis (H): \"{hypothesis}\"")
    print("-" * 40)

    print("Step 1: Find the relation between the Premise (P) and the un-negated Hypothesis (H').")
    print(f"H' = \"{unnegated_hypothesis}\"")
    print("\nFirst, we find the relation for the edits transforming P to H':")
    
    # Edit 1: Deletion
    edit1_relation = "<"
    edit1_name = "Forward Entailment"
    print(f"  - Edit 'pop song' -> 'song' has the relation: {edit1_relation} ({edit1_name})")
    
    # Edit 2: Substitution
    edit2_relation = "|"
    edit2_name = "Alternation"
    print(f"  - Edit 'Taylor Swift' -> 'Michael Jackson' has the relation: {edit2_relation} ({edit2_name})")

    # Composition for H'
    rel_p_h_prime = "|"
    rel_p_h_prime_name = "Alternation"
    print(f"\nComposing these edits using the join table (join({edit1_relation}, {edit2_relation})) gives the relation between P and H':")
    print(f"Rel(P, H') = {rel_p_h_prime} ({rel_p_h_prime_name})")
    print("-" * 40)

    print("Step 2: Compose the result from Step 1 with the negation.")
    negation_relation = "^"
    negation_name = "Negation"
    print(f"The relation between H' and H (which is 'not H') is {negation_relation} ({negation_name}).")
    
    # Final Composition
    final_relation = "<"
    final_relation_name = "Forward Entailment"
    print("\nUsing the composition table for the final step:")
    print(f"Final Relation = compose(Rel(P, H'), {negation_relation})")
    print(f"Final Equation: compose({rel_p_h_prime}, {negation_relation}) = {final_relation}")
    print("-" * 40)
    
    print(f"The final projected natural logic operator is {final_relation}, which is named {final_relation_name}.")

solve_entailment()