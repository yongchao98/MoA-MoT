def solve_entailment():
    """
    Determines the natural logic operator for the given premise and hypothesis
    by modeling the inference as a chain of atomic steps and composing their relations.
    """

    # Define the propositions in the logical chain from Premise to Hypothesis.
    premise = "Mark is singing a pop song by Taylor Swift"
    p_intermediate_1 = "Mark is singing a song by Taylor Swift"
    p_intermediate_2 = "Mark is singing a song by Michael Jackson"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    # Step 1: Define the relation between the Premise and the first intermediate proposition.
    # A "pop song" is a type of "song". If the premise is true, p_intermediate_1 must be true.
    # This is Forward Entailment (<).
    rel_premise_to_p1 = "Forward Entailment"

    # Step 2: Define the relation between the first and second intermediate propositions.
    # "Taylor Swift" and "Michael Jackson" are disjoint artists. A single song cannot be by both.
    # Therefore, the two propositions are mutually exclusive.
    # This is Alternation (|).
    rel_p1_to_p2 = "Alternation"

    # Step 3: Define the relation between the second intermediate proposition and the Hypothesis.
    # The hypothesis is the direct negation of the second intermediate proposition.
    # This is Negation (^).
    rel_p2_to_hypothesis = "Negation"

    # Now, we compose these relations to find the final relation between the Premise and Hypothesis.
    # We'll define a composition function based on logical deduction.
    # e.g., compose(rel(A, B), rel(B, C)) -> rel(A, C)
    def compose(rel1, rel2):
        """Composes two sequential logical relations."""
        # Composition rule 1: Forward Entailment (<) followed by Alternation (|)
        # If A entails B (A => B), and B and C are alternates (B => not C),
        # then by chaining the logic, A entails not C (A => not C).
        # The relationship between A and C is therefore Negation (^).
        if rel1 == "Forward Entailment" and rel2 == "Alternation":
            return "Negation"

        # Composition rule 2: Negation (^) followed by Negation (^)
        # If A and B are contradictory (A => not B), and B and C are contradictory (B <=> not C),
        # then by substitution, we get A => not(not C), which simplifies to A => C.
        # The relationship between A and C is therefore Forward Entailment (<).
        if rel1 == "Negation" and rel2 == "Negation":
            return "Forward Entailment"
        
        return "Unknown Composition"

    # First composition: Find the relation between the Premise and p_intermediate_2.
    rel_premise_to_p2 = compose(rel_premise_to_p1, rel_p1_to_p2)

    # Final composition: Find the relation between the Premise and the Hypothesis.
    final_relation = compose(rel_premise_to_p2, rel_p2_to_hypothesis)

    # --- Output the reasoning and the final result ---
    print("Step-by-step derivation of the entailment relation:")
    print("--------------------------------------------------")
    print(f"Premise (P1): '{premise}'")
    print(f"Hypothesis (P4): '{hypothesis}'\n")

    print("The inference is modeled as a chain of propositions:")
    print(f"P1 -> P2: From '{premise}' to '{p_intermediate_1}'")
    print(f"   - Relation: {rel_premise_to_p1} (<)\n")

    print(f"P2 -> P3: From '{p_intermediate_1}' to '{p_intermediate_2}'")
    print(f"   - Relation: {rel_p1_to_p2} (|)\n")
    
    print(f"P3 -> P4: From '{p_intermediate_2}' to '{hypothesis}'")
    print(f"   - Relation: {rel_p2_to_hypothesis} (^)\n")

    print("Composing the relations to find the final relation between Premise (P1) and Hypothesis (P4):")
    print("-----------------------------------------------------------------------------------------")
    print(f"Composition 1: Relation(P1, P3) = Compose(Relation(P1, P2), Relation(P2, P3))")
    print(f"   - Compose('{rel_premise_to_p1}', '{rel_p1_to_p2}') => '{rel_premise_to_p2}'\n")

    print(f"Composition 2: Relation(P1, P4) = Compose(Relation(P1, P3), Relation(P3, P4))")
    print(f"   - Compose('{rel_premise_to_p2}', '{rel_p2_to_hypothesis}') => '{final_relation}'\n")

    print(f"The final projected natural logic operator is: {final_relation}")

solve_entailment()