def check_betti_number_claim(n):
    """
    This function checks the claim in question (b) for a given n.

    The claim is that for G = SU(n), the second Betti number b2 of any
    coadjoint orbit is always n - 1.

    We test this by comparing the b2 for a generic orbit (the full flag
    manifold) with the b2 for a specific singular orbit (a Grassmannian).
    """

    print(f"--- Checking claim for G = SU({n}) ---")

    # For a generic (regular) coadjoint orbit, O_lambda = SU(n)/T,
    # where T is the maximal torus.
    # The rank of SU(n) is n-1.
    # The second Betti number b2(SU(n)/T) is equal to the rank of SU(n).
    b2_generic = n - 1
    print(f"For a generic orbit, the expected second Betti number is n - 1.")
    print(f"Equation: b2_generic = n - 1 = {n} - 1 = {b2_generic}")

    # For a singular coadjoint orbit, we can take the Grassmannian Gr(k, n).
    # Let's choose k=2 (for n > 2).
    # Gr(k, C^n) is the space of k-dimensional complex subspaces in C^n.
    # It is a well-known result that its second Betti number is 1.
    # This orbit corresponds to a specific choice of singular lambda.
    if n > 2:
        k = 2
        b2_singular = 1
        print(f"\nNow consider a singular orbit, the Grassmannian Gr({k}, {n}).")
        print(f"The second Betti number for this orbit is known to be 1.")
        print(f"Equation: b2_singular = 1")

        # Compare the results
        print("\n--- Conclusion ---")
        print(f"Is b2 always equal to n - 1?")
        is_always_true = (b2_generic == b2_singular)
        if not is_always_true:
            print(f"No. For G = SU({n}), the value for a generic orbit is {b2_generic}, but for the singular orbit Gr({k},{n}) it is {b2_singular}.")
            print(f"The statement is false because {b2_singular} != {b2_generic}.")
        else:
            # This would only happen for n=2, where Gr(1,2)=CP^1 and n-1=1.
            print(f"Yes, for n={n}, both values are {b2_generic}.")
            print("However, the claim must hold for ALL n. Let's try a larger n.")

# Let's run the check for n=4, which provides a clear counterexample.
check_betti_number_claim(4)
