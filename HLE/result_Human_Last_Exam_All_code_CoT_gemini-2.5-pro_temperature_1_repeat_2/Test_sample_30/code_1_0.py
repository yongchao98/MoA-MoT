def explain_correct_statement():
    """
    This function prints the correct statement about interactive proof systems and explains why it is correct.
    """
    correct_option = "D"
    statement_text = "If a prover and verifier are required, there are exponential requirements on the computational power of the prover, whereas the verifier is required to run in polynomial time."

    explanation = """
This statement accurately describes the fundamental model of an interactive proof system.
1.  **Efficient Verifier**: The 'proof' is only useful if it can be checked quickly. In computational complexity, 'quickly' or 'efficiently' is formalized as running in polynomial time with respect to the input size. This ensures the verifier is a practical entity.
2.  **Powerful Prover**: The prover's task is to provide the information that convinces the verifier. For complex problems (e.g., those complete for PSPACE), finding this information can be computationally very hard. Therefore, the prover is modeled as having immense computational resources, often super-polynomial (like exponential time) or even being computationally unbounded.

This contrast between a computationally limited verifier and a powerful (but untrusted) prover is what gives interactive proof systems their unique power, allowing them to solve problems beyond the scope of static proofs (NP).
"""

    print(f"The correct statement is option {correct_option}:")
    print(f'"{statement_text}"')
    print("\n--- Explanation ---")
    print(explanation)

explain_correct_statement()