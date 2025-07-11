def solve_natural_logic_inference():
    """
    Calculates and explains the final natural logic operator for a given premise and hypothesis
    based on MacCartney's compositional framework.
    """
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print(f"Premise (P): \"{premise}\"")
    print(f"Hypothesis (H): \"{hypothesis}\"")
    print("\nWe will find the relation R(P, H) by composing the relations of the edits that transform P into H.\n")

    # Step-by-step derivation
    print("--- Derivation Steps ---")

    # Step 0
    print("Step 0: Initial State")
    print("We start with the Premise compared to itself.")
    print("Relation: P vs. P")
    print("Operator: EQUIVALENCE (=)\n")

    # Step 1
    print("Step 1: Insertion of 'not'")
    print("The first edit changes 'is singing' to 'is NOT singing'.")
    print("The semantic relation for this negation edit is ALTERNATION (|).")
    print("We compose this with our starting relation:")
    print("Intermediate Relation = Join(EQUIVALENCE, ALTERNATION) = ALTERNATION (|)\n")
    
    # Step 2
    print("Step 2: Deletion of 'pop'")
    print("The second edit changes '...not singing a pop song...' to '...not singing a song...'.")
    print("The relation of this transition within a negative context is REVERSE ENTAILMENT (^).")
    print("   (Because 'not singing a song' implies 'not singing a pop song').")
    print("We compose this with our current relation:")
    print("Intermediate Relation = Join(ALTERNATION, REVERSE ENTAILMENT) = CONTRADICTION (~)\n")

    # Step 3
    print("Step 3: Substitution of 'Taylor Swift' with 'Michael Jackson'")
    print("The final edit changes '...song by Taylor Swift' to '...song by Michael Jackson'.")
    print("The relation of this transition ('...not...by TS' -> '...not...by MJ') is COVER (v).")
    print("   (Because if he IS singing a song by TS, he is NOT singing a song by MJ).")
    print("We compose this with our current relation to get the final result:")
    print("Final Relation = Join(CONTRADICTION, COVER) = FORWARD ENTAILMENT (ㄷ)\n")
    
    # Final Answer
    print("--- Conclusion ---")
    final_operator = "FORWARD ENTAILMENT"
    print("The final projected natural logic operator is FORWARD ENTAILMENT.")
    
    # As requested, printing the elements of the final equation
    print("\nFinal compositional step:")
    part1 = "CONTRADICTION"
    part2 = "COVER"
    result = "FORWARD ENTAILMENT"
    print(f"Join({part1}, {part2}) = {result}")

solve_natural_logic_inference()