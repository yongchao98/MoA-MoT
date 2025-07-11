def solve_set_theory_question():
    """
    Analyzes the logical relationship between large cardinal axioms to answer the question.
    """
    print("Step 1: Define the components of the problem.")
    system_S = "ZFC + 'There exists an inaccessible cardinal κ'"
    statement_P_definition = "'There exists a nontrivial elementary embedding j: V → M with critical point κ'"
    statement_P_equivalence = "'κ is a measurable cardinal'"
    print(f"System S is: {system_S}")
    print(f"Statement P is: {statement_P_definition}, which is equivalent to {statement_P_equivalence}.")
    print("-" * 30)

    print("Step 2: Compare the consistency strength of the axioms.")
    print("In set theory, the existence of a measurable cardinal is a strictly stronger axiom than the existence of an inaccessible cardinal.")
    print("A measurable cardinal is provably inaccessible, but an inaccessible cardinal is NOT provably measurable.")
    print("-" * 30)

    print("Step 3: Analyze the provability of P in S.")
    print("Could S prove P? No. Proving P from S would mean proving a stronger theory from a weaker one, which is not possible.")
    print("This eliminates choice A.")
    print("-" * 30)

    print("Step 4: Analyze the disprovability of P in S.")
    print("Could S prove not(P)? No. Proving not(P) would mean the existence of an inaccessible cardinal contradicts the existence of a measurable cardinal.")
    print("This is false, as the theory ZFC + 'a measurable cardinal exists' is considered consistent.")
    print("This eliminates choices B and D.")
    print("-" * 30)

    print("Step 5: Conclude based on the analysis.")
    print("Since P can neither be proved nor disproved from S, P is independent of S.")
    print("The correct choice is C.")
    print("-" * 30)

    # Final Answer
    final_answer = 'C'
    print(f"The final answer is: {final_answer}")


solve_set_theory_question()