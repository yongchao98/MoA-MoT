def check_betti_number_claim():
    """
    This function provides a specific counterexample to the claim in part (b).

    The claim is: For G = SU(n), with lambda in the Weyl alcove,
    is the second Betti number b_2(O_lambda) always given by n - 1?
    """

    # We select n=3 as our case study.
    n = 3
    # The rank of SU(n) is n-1.
    rank = n - 1
    
    print("--- Analyzing the claim for G = SU(n) ---")
    print(f"Let's test the case where n = {n}.")
    print(f"The rank of SU({n}) is n - 1 = {rank}.")
    print(f"The claim is that for any coadjoint orbit O_lambda of SU({n}), the second Betti number is b_2 = {rank}.")
    print("-" * 20)

    # We consider a coadjoint orbit corresponding to a "singular" element lambda.
    # This orbit is the complex projective plane, CP^2.
    orbit_name = "CP^2"
    
    # The Betti numbers for CP^m are well-established. b_2(CP^m) is always 1 for m>=1.
    # For CP^2 (where m=2), b_2 is 1.
    b2_of_orbit = 1

    print(f"We consider a specific coadjoint orbit O_lambda, which is the complex projective plane, {orbit_name}.")
    print(f"The second Betti number for {orbit_name} is known to be b_2 = {b2_of_orbit}.")
    
    # We now compare the actual b_2 with the claimed value.
    print(f"\nComparing the orbit's b_2 with the group's rank:")
    print(f"Actual b_2 for {orbit_name}: {b2_of_orbit}")
    print(f"Claimed value (rank of SU({n})): {rank}")

    if b2_of_orbit == rank:
        result = "The claim holds for this example."
        final_answer = "Yes"
    else:
        result = "The claim does not hold for this example."
        final_answer = "No"

    print(f"\nResult: {result}")
    print("\n--- Conclusion for Part (b) ---")
    print(f"We have found a coadjoint orbit of SU({n}) whose second Betti number ({b2_of_orbit})")
    print(f"is not equal to the rank of the group ({rank}).")
    print(f"Therefore, the statement 'b_2(O_lambda) is always n - 1' is false.")
    print(f"Final Answer for (b): {final_answer}")


if __name__ == '__main__':
    check_betti_number_claim()
