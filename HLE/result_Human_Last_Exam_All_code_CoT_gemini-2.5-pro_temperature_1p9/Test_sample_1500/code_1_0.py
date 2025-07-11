def check_betti_number_condition(n, k):
    """
    This function demonstrates the counterexample for question (b).
    It checks if the second Betti number of the Grassmannian Gr(k,n),
    which is a coadjoint orbit of SU(n), is equal to n-1.
    """
    # The complex Grassmannian Gr(k,n) is a coadjoint orbit of SU(n).
    # For 0 < k < n, its second Betti number b_2 is known to be 1.
    # The rank of SU(n) is n-1. The question is if b_2 is always n-1.

    if not (isinstance(n, int) and isinstance(k, int) and 0 < k < n):
        print(f"Error: Please provide integers n, k such that 0 < k < n.")
        return

    b2_Gr = 1
    rank_SU_n = n - 1

    print(f"--- Counterexample for part (b) using G = SU({n}) ---")
    print(f"Consider the coadjoint orbit corresponding to the Grassmannian Gr({k}, {n}).")
    print(f"The second Betti number of this orbit is b_2(Gr({k}, {n})) = {b2_Gr}.")
    print(f"The rank of SU({n}), which is n-1, is {rank_SU_n}.")
    
    is_equal = (b2_Gr == rank_SU_n)
    print(f"Is b_2 equal to n-1? {is_equal}")

    if not is_equal:
        print(f"Result: {b2_Gr} is not equal to {rank_SU_n}.")
        print("This confirms that b_2 is not always equal to n-1 for all coadjoint orbits.")
    else:
        print(f"Result: {b2_Gr} is equal to {rank_SU_n}.")
        print("This case does not serve as a counterexample. Try n>2.")

# We run the check for SU(4) and the orbit Gr(2, 4).
check_betti_number_condition(4, 2)