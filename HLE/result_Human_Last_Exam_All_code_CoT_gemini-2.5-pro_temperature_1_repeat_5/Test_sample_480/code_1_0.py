def solve_natural_logic_inference():
    """
    Calculates the final projected natural logic operator for the given
    premise and hypothesis by composing relations from a series of edits.
    """

    # Premise (P): "Mark is singing a pop song by Taylor Swift"
    # Hypothesis (H): "Mark is not singing a song by Michael Jackson"
    # Path: P -> H1 -> H2 -> H

    print("Deriving the final relation step-by-step:")
    print("=" * 40)

    # --- Step 1: P -> H1 ---
    # Edit: Delete "pop". P becomes H1: "Mark is singing a song by Taylor Swift"
    # Lexical relation ('pop song' vs 'song') is Forward Entailment (f).
    # The context is positive (upward monotonic), so the projected relation is f.
    relation_step1 = "f"
    print("Step 1: Delete 'pop'")
    print(f"Relation R(P, H1): {relation_step1} (Forward Entailment)")
    print("-" * 40)

    # --- Step 2: H1 -> H2 ---
    # Edit: Substitute "Taylor Swift" with "Michael Jackson".
    # H1 becomes H2: "Mark is singing a song by Michael Jackson"
    # Lexical relation ('Taylor Swift' vs 'Michael Jackson') is Exclusion (e).
    # The context is positive (upward monotonic), so the projected relation is e.
    relation_step2 = "e"
    print("Step 2: 'Taylor Swift' -> 'Michael Jackson'")
    print(f"Relation R(H1, H2): {relation_step2} (Exclusion)")
    print("-" * 40)

    # --- Step 3: H2 -> H ---
    # Edit: Negate the sentence.
    # H2 becomes H: "Mark is not singing a song by Michael Jackson"
    # The relation between a sentence and its negation is Alternation (a).
    relation_step3 = "a"
    print("Step 3: Negate the sentence")
    print(f"Relation R(H2, H): {relation_step3} (Alternation)")
    print("-" * 40)

    # --- Composition ---
    # We compose the relations: R(P,H) = R(P,H1) ; R(H1,H2) ; R(H2,H)

    # First composition: R(P, H2) = R(P,H1) ; R(H1,H2)
    # According to MacCartney's table, f ; e = e
    intermediate_relation = "e"
    print("Composition 1: R(P, H1) ; R(H1, H2)")
    print(f"Equation: {relation_step1} ; {relation_step2} = {intermediate_relation}")
    print(f"Intermediate Relation R(P, H2) is: {intermediate_relation} (Exclusion)")
    print("-" * 40)

    # Final composition: R(P, H) = R(P,H2) ; R(H2,H)
    # According to MacCartney's table, e ; a = f
    final_relation = "f"
    print("Composition 2 (Final): R(P, H2) ; R(H2, H)")
    print(f"Equation: {intermediate_relation} ; {relation_step3} = {final_relation}")
    print("-" * 40)

    final_operator_name = "Forward Entailment"
    print(f"The final projected natural logic operator is '{final_relation}'.")
    print(f"The name of the operator is: {final_operator_name}")

solve_natural_logic_inference()