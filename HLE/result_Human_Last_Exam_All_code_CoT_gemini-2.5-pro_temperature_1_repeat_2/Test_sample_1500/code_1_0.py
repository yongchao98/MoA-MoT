def check_betti_number_for_su_n(n):
    """
    Checks if the second Betti number b2 of a coadjoint orbit of SU(n)
    is always equal to n-1 by comparing a regular and a singular orbit.

    Args:
        n (int): The dimension for the special unitary group SU(n).
                 Should be an integer greater than 1.
    """
    if not isinstance(n, int) or n < 2:
        print("Please provide an integer n >= 2.")
        return

    # Case 1: Lambda is regular (in the interior of the Weyl alcove).
    # The orbit is the full flag manifold SU(n)/T.
    # The second Betti number b2 is the rank of SU(n), which is n-1.
    b2_regular = n - 1

    # Case 2: Lambda is singular (on the boundary of the Weyl alcove).
    # A possible orbit is the complex projective space CP^{n-1}.
    # The second Betti number b2 of CP^{m} for m>=1 is 1.
    b2_singular = 1

    print(f"Analysis for G = SU({n}):")
    print("=" * 30)
    print("Question: Is the second Betti number b2(O_lambda) always n-1?")
    print("")

    print("Consider a regular orbit (lambda in interior of Weyl alcove):")
    print(f"  - Orbit type: Full flag manifold SU({n})/T")
    print(f"  - The second Betti number is the rank of SU({n}).")
    print(f"  - b2 = n - 1")
    print(f"  - b2 = {n} - 1 = {b2_regular}")
    print("")

    print("Consider a singular orbit (lambda on boundary of Weyl alcove):")
    print(f"  - Orbit type: Complex projective space CP^{n-1} = CP^{{{n-1}}}")
    print(f"  - The second Betti number of CP^{{{n-1}}} is 1.")
    print(f"  - b2 = {b2_singular}")
    print("")

    print("Conclusion:")
    if b2_regular == b2_singular:
        print(f"For n={n}, the values are equal ({b2_regular} == {b2_singular}).")
        print("The statement holds in this case.")
    else:
        print(f"For n={n}, the values are NOT equal ({b2_regular} != {b2_singular}).")
        print("This provides a counterexample to the claim that b2 is *always* n-1.")
        print("Therefore, the answer to question (b) is 'No'.")

# Run the check for n=4 as a demonstration.
check_betti_number_for_su_n(4)
