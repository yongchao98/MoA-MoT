def verify_betti_number_claim():
    """
    This script verifies the claim in part (b) for G = SU(n).
    The claim is that the second Betti number b_2(O_lambda) is always n-1.
    We will test this claim with a counterexample for n=3.
    """
    n = 3
    n_minus_1 = n - 1

    print(f"We are investigating the claim for G = SU(n) with n = {n}.")
    print(f"The claim states that the second Betti number b_2 is always equal to n-1, which is {n_minus_1}.")
    print("-" * 50)

    # The rank of SU(n) is n-1. For SU(3), the rank is 2.
    # The number of simple roots of SU(3) is 2.
    rank_G = n - 1
    
    # --- Case 1: lambda is regular ---
    # In this case, the orbit is the full flag manifold SU(3)/T.
    # The stabilizer is the maximal torus T, whose semisimple part has rank 0.
    rank_stabilizer_regular = 0
    b2_regular = rank_G - rank_stabilizer_regular
    
    print("Case 1: lambda is regular (in the interior of the Weyl alcove)")
    print("The stabilizer is the maximal torus T.")
    print(f"The number of simple roots of the stabilizer's semisimple part is {rank_stabilizer_regular}.")
    print("The second Betti number b_2 is the number of simple roots of the group minus those of the stabilizer.")
    print(f"b_2 = {rank_G} - {rank_stabilizer_regular} = {b2_regular}")
    print(f"Here, b_2 = {b2_regular} is equal to n-1 = {n_minus_1}. The claim holds.")
    print("-" * 50)

    # --- Case 2: lambda is singular ---
    # For a singular lambda, the stabilizer is a larger Levi subgroup.
    # For SU(3), a singular lambda can have a stabilizer L = S(U(2) x U(1)).
    # The semisimple part of L is SU(2), which has rank 1.
    rank_stabilizer_singular = 1
    b2_singular = rank_G - rank_stabilizer_singular

    print("Case 2: lambda is singular (on a wall of the Weyl alcove)")
    print("A possible stabilizer is L = S(U(2) x U(1)).")
    print(f"The semisimple part of this stabilizer is SU(2), which has a root system of rank {rank_stabilizer_singular}.")
    print("The second Betti number b_2 is calculated similarly.")
    print(f"b_2 = {rank_G} - {rank_stabilizer_singular} = {b2_singular}")
    print(f"Here, b_2 = {b2_singular} is NOT equal to n-1 = {n_minus_1}. The claim fails.")
    print("-" * 50)

    print("Conclusion: As we have found a counterexample, the statement that b_2(O_lambda) is *always* n-1 is false.")

verify_betti_number_claim()