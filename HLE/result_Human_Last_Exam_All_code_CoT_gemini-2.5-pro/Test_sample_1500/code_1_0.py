def check_betti_number_statement(n):
    """
    This function tests the statement from question (b) for G = SU(n).
    The statement claims that the second Betti number b_2 is always n - 1.
    We use a known counterexample: the coadjoint orbit O_omega_1, which is
    isomorphic to the complex projective space CP^(n-1). The second Betti
    number of CP^(n-1) is 1 for all n >= 2.
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    # The actual second Betti number for the orbit O_omega_1 ~ CP^(n-1)
    b2_orbit = 1
    
    # The value predicted by the statement in the question
    predicted_b2 = n - 1
    
    print(f"--- Analysis for G = SU({n}) ---")
    print(f"Consider the coadjoint orbit corresponding to the fundamental weight omega_1.")
    print(f"This orbit is the complex projective space CP^({n-1}).")
    print(f"The actual second Betti number is b_2(CP^{n-1}) = {b2_orbit}.")
    print(f"The statement in question (b) suggests the value should be n - 1.")
    print(f"For n = {n}, the suggested value is {n} - 1 = {predicted_b2}.")
    
    # Final comparison
    print(f"Comparing the actual value to the suggested value: {b2_orbit} == {predicted_b2} is {b2_orbit == predicted_b2}.")
    if b2_orbit != predicted_b2:
        print("Since the values do not match, the statement is false.")

# We demonstrate the counterexample with n=3.
check_betti_number_statement(3)