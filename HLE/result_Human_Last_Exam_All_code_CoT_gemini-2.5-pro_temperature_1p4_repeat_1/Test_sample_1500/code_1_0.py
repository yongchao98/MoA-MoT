def check_b2_for_sun(n):
    """
    Checks if the second Betti number b_2 for coadjoint orbits of SU(n)
    is always equal to n-1. This script provides a counterexample for n > 2.

    Args:
        n (int): The dimension for the special unitary group SU(n).
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    rank = n - 1
    print(f"Let's analyze the case for G = SU(n) where n = {n}.")
    print(f"The rank of SU({n}) is n-1, which is {rank}.")
    print("The question is whether the second Betti number b_2 of any coadjoint orbit is always equal to this rank.")
    print("-" * 30)

    # Case 1: Regular orbit
    # For a regular element lambda, the orbit is the full flag manifold SU(n)/T.
    # Its second Betti number b_2 is equal to the rank of the group.
    b2_regular = rank
    print("Consider a 'regular' orbit, where the stabilizer is the maximal torus T.")
    print("This orbit is the full flag manifold SU(n)/T.")
    print(f"For this orbit, the second Betti number b_2 is equal to the rank of SU({n}).")
    print(f"So, b_2 = {n} - 1 = {b2_regular}.")
    print("-" * 30)

    # Case 2: A specific singular orbit
    # Consider the orbit that is diffeomorphic to the complex projective space CP^{n-1}.
    # The second Betti number of CP^{k} is 1 for k >= 1.
    b2_singular = 1
    print("Now, consider a 'singular' orbit, specifically the one diffeomorphic to the complex projective space CP^{n-1}.")
    print(f"For n={n}, this orbit is the complex projective space CP^{n-1} = CP^{n - 1} = CP^{2}.")
    print("The second Betti number of any complex projective space CP^k (for k>=1) is 1.")
    print(f"So, for this orbit, b_2 = {b2_singular}.")
    print("-" * 30)

    # Compare the results and conclude
    print("Conclusion for question (b):")
    print("We have found two different coadjoint orbits for SU(n) with different second Betti numbers:")
    print(f"  - For the regular orbit, b_2 = {b2_regular}")
    print(f"  - For the singular orbit CP^{n-1}, b_2 = {b2_singular}")
    print(f"Is {b2_regular} always equal to {b2_singular}?")
    if b2_regular == b2_singular:
        print(f"For n={n}, the values happen to be equal ({b2_regular} == {b2_singular}). This occurs only when n=2.")
    else:
        print(f"For n={n}, the values are not equal: {b2_regular} != {b2_singular}.")
    
    print("\nThis provides a clear counterexample. Therefore, the statement that b_2(O_lambda) is *always* n-1 is false.")

# Run the demonstration for n=3, a case where the values are different.
check_b2_for_sun(3)