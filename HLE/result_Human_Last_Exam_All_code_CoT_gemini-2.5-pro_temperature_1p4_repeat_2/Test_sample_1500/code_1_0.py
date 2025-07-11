def check_betti_number_statement(n):
    """
    This function provides a counterexample to the statement in question (b).

    It considers the coadjoint orbit of SU(n) diffeomorphic to the complex
    projective space CP^{n-1}. It then compares the actual second Betti number
    of this space with the proposed formula n-1.
    """
    if not isinstance(n, int) or n < 2:
        print("Please provide an integer n >= 2.")
        return

    # The complex projective space CP^{n-1} is a coadjoint orbit for SU(n).
    # Its Betti numbers are b_0=1, b_2=1, ..., b_{2(n-1)}=1, and others 0.
    # Therefore, the second Betti number is always 1.
    b2_orbit = 1

    # The value proposed in the question.
    formula_val = n - 1

    print(f"--- Checking statement (b) for G = SU({n}) ---")
    print(f"Consider the coadjoint orbit O_lambda which is the complex projective space CP^{n-1}.")
    print(f"The second Betti number for this orbit is b_2(CP^{n-1}) = {b2_orbit}.")
    
    # The final equation mentioned in the problem description is the comparison.
    # We output each number in this final comparison.
    print(f"The value from the statement is n - 1 = {n} - 1 = {formula_val}.")

    if b2_orbit == formula_val:
        print(f"For n={n}, the values are equal: {b2_orbit} == {formula_val}.")
        print("However, the statement must hold for ALL orbits and all n.")
    else:
        print(f"For n={n}, the values are NOT equal: {b2_orbit} != {formula_val}.")
        print("This serves as a counterexample, so statement (b) is false.")
    
    print("-" * 40)


# Run the check for a specific case, e.g., n=3.
check_betti_number_statement(3)

# The final answers based on the reasoning provided above.
final_answer = "(a) True; (b) No; (c) No"
print("\nFinal Answer:")
print(final_answer)
print("<<<" + final_answer + ">>>")
