def check_betti_number_for_su_n():
    """
    This function provides a specific counterexample to the statement in question (b).
    The question is: For G = SU(n), is the second Betti number b_2(O_lambda) always n - 1?
    We will test this for n=3.
    """

    n = 3
    group = "SU(3)"

    # According to the statement, the Betti number should be n-1.
    proposed_b2 = n - 1
    
    # However, for G = SU(3), there exists a coadjoint orbit O_lambda
    # that is diffeomorphic to the complex projective space CP^2.
    manifold_example = "CP^2"
    
    # The second Betti number of CP^k is 1 for any k. So for CP^2, it is 1.
    actual_b2 = 1

    print(f"Investigating the claim for G = {group} (n={n}).")
    print(f"The proposed formula states that the second Betti number should be b_2 = n - 1.")
    print(f"For n = {n}, this gives the value: {proposed_b2}")
    
    print(f"\nLet's consider a counterexample.")
    print(f"A coadjoint orbit for {group} is the complex projective space {manifold_example}.")
    print(f"The second Betti number of {manifold_example} is b_2({manifold_example}) = {actual_b2}.")

    print(f"\nNow we compare the proposed value with the actual value for this example.")
    print(f"Is {actual_b2} == {n} - 1?")
    print(f"Final Equation: {actual_b2} == {proposed_b2}")
    print(f"Result: {actual_b2 == proposed_b2}")
    print("\nSince the values are not equal, the statement is false.")

check_betti_number_for_su_n()