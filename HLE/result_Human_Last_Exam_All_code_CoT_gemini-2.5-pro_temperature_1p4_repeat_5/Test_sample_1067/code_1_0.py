def solve_logic_puzzle():
    """
    This function explains the logical deduction for the chosen answer.
    It demonstrates why option C provides a valid explanation for the apparent contradiction.
    """
    
    # Define the propositions from Choice C
    P = "The dog detects an intruder"
    Q = "The dog barked"
    R = "The dog was asleep"

    # State the known facts from the problem description
    fact_P = True  # "The dog detected an intruder"
    fact_not_Q = True  # "The dog did not bark"

    print("Analyzing the Logic of Choice C:")
    print("-" * 30)
    print(f"P: {P}")
    print(f"Q: {Q}")
    print(f"R: {R}\n")

    print("The argument presented is: [(P ∧ ¬R) → Q] ∧ (¬Q ∧ P), ∴ R\n")

    print("Step-by-step deduction:")
    print("1. We accept Premise 1, a refined rule: (P ∧ ¬R) → Q")
    print("   Which means: 'If the dog detects an intruder AND is NOT asleep, then it will bark.'")
    
    print("\n2. We are given Premise 2 from the verifiable proof: ¬Q ∧ P")
    print("   Which means: 'The dog did not bark' (¬Q is True) AND 'The dog detected an intruder' (P is True).")

    print("\n3. Let's use the contrapositive of Premise 1: ¬Q → ¬(P ∧ ¬R)")
    print("   Which means: 'If the dog did not bark, then it is NOT the case that (the dog detected an intruder AND was not asleep)'.")

    print("\n4. Since '¬Q' is True (from Step 2), we can conclude by Modus Ponens that ¬(P ∧ ¬R) must be True.")

    print("\n5. Applying De Morgan's laws to ¬(P ∧ ¬R), we get (¬P ∨ R).")
    print("   Which means: 'Either the dog did not detect an intruder, OR the dog was asleep'.")

    print("\n6. From Step 2, we know that P is True, which means ¬P must be False.")
    
    print("\n7. In the statement (¬P ∨ R), since ¬P is False, R must be True for the whole statement to be True.")
    print("   This form of reasoning is called a Disjunctive Syllogism.")

    print("\nConclusion: The deduction is logically valid. The reason the dog did not bark is R: The dog was asleep.")
    print("-" * 30)

    # Final Answer Formatting
    print("The final logical statement is:")
    p_statement = "P: The dog detects an intruder"
    q_statement = "Q: The dog barked"
    r_statement = "R: The dog was asleep"
    
    final_equation = "[(P ∧ ¬R)→Q] ∧ (¬Q ∧ P),∴R"

    print(f"{p_statement}")
    print(f"{q_statement}")
    print(f"{r_statement}")
    print(f"\n{final_equation}")


solve_logic_puzzle()