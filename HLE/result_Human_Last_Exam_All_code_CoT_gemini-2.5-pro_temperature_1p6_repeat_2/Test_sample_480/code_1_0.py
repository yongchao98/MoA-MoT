def solve_natural_logic_inference():
    """
    This function presents the step-by-step derivation of the final natural
    logic operator for the given premise and hypothesis, according to
    MacCartney's compositional framework.
    """
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    # The following variables store the results of the logical analysis explained above.
    
    # Starting point: Relation between a sentence and itself is Identity (=).
    # Edit 1 (insert 'not'): Changes relation to Negation (^).
    # Composition 1: T(=, ^) = ^
    r1 = "^"
    r1_name = "Negation"

    # Edit 2 (delete 'pop' in downward context): Projected relation is < (new sentence is entailed by old).
    # Composition 2: T(^, <) = |
    r2 = "|"
    r2_name = "Alternation"

    # Edit 3 (substitute 'TS' with 'MJ' in downward context): Projected relation is |.
    # Composition 3: T(|, |) = |
    final_op_symbol = "|"
    final_op_name = "Alternation"

    # --- Outputting the results ---
    print("Deriving the final projected natural logic operator:")
    print("-" * 50)
    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"")
    print("-" * 50)

    print("The derivation involves composing the semantic relations of three left-to-right edits:")
    
    print("\n1. Insertion of 'not':")
    print(f"   - The initial relation becomes {r1} ({r1_name}).")

    print("\n2. Deletion of 'pop' in a negative context:")
    print("   - The composition T(Negation, Forward Entailment) is performed.")
    print(f"   - The intermediate relation becomes {r2} ({r2_name}).")

    print("\n3. Substitution of 'Taylor Swift' with 'Michael Jackson' in a negative context:")
    print(f"   - The final composition T(Alternation, Alternation) is performed.")
    print(f"   - The final relation is determined to be {final_op_symbol} ({final_op_name}).")
    
    print("\n" + "-" * 50)
    print("The final equation showing the relationship is:")
    print(f"\"{premise}\"")
    print(f"    {final_op_symbol}")
    print(f"\"{hypothesis}\"")
    
    print(f"\nThe name of the final projected natural logic operator is {final_op_name}.")


solve_natural_logic_inference()