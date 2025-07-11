#
# This script verifies the answer to part (b) for G = SU(n).
# It shows that the second Betti number b2 is not always n-1.
#

def verify_b2_for_sun(n):
    """
    Checks if b2 is always n-1 for coadjoint orbits of SU(n).
    
    Args:
        n (int): The dimension for the SU(n) group. Must be >= 3 for a good counterexample.
    """
    if n < 3:
        print(f"For SU({n}), the non-regular cases are trivial. Please use n >= 3.")
        return

    print(f"--- Verification for G = SU({n}) ---")
    print(f"The question is whether the second Betti number b2 of any coadjoint orbit is always equal to n - 1 = {n - 1}.")
    print("\nCoadjoint orbits of SU(n) correspond to flag manifolds F(p_1, ..., p_k),")
    print("where (p_1, ..., p_k) is an integer partition of n.")
    print("The second Betti number, b2, for such a manifold is given by the formula: b2 = k - 1.")

    # Case 1: The regular orbit
    # This corresponds to the partition where every part is 1.
    partition_regular = tuple([1] * n)
    k_regular = len(partition_regular)
    b2_regular = k_regular - 1
    
    print(f"\n1. Consider a REGULAR orbit.")
    print(f"This corresponds to the partition {partition_regular}, which gives the full flag manifold SU({n})/T.")
    print(f"For this partition, the number of parts is k = {k_regular}.")
    print(f"The second Betti number is b2 = k - 1 = {k_regular} - 1 = {b2_regular}.")
    print(f"The proposed value from the question is n - 1 = {n} - 1 = {n-1}.")
    print(f"In this case, the values match: {b2_regular} == {n-1}. The statement holds for regular orbits.")

    # Case 2: A singular orbit (provides a counterexample)
    # For n=3, a simple singular partition is (2, 1).
    partition_singular = (n - 1, 1)
    k_singular = len(partition_singular)
    b2_singular = k_singular - 1

    print(f"\n2. Now consider a SINGULAR orbit.")
    print(f"An example is the orbit corresponding to the partition {partition_singular}.")
    print(f"This orbit is the projective space CP^{n-1}.")
    print(f"For this partition, the number of parts is k = {k_singular}.")
    print(f"The second Betti number is b2 = k - 1 = {k_singular} - 1 = {b2_singular}.")
    print(f"The proposed value is still n - 1 = {n} - 1 = {n-1}.")

    # Final conclusion
    print("\n--- Conclusion ---")
    print("Comparing the calculated b2 for the singular orbit with the proposed value:")
    print(f"Calculated b2 = {b2_singular}")
    print(f"Proposed value n-1 = {n-1}")
    print(f"Since {b2_singular} != {n-1}, the statement that b2 is ALWAYS n-1 is false.")

# Run the verification for n=3 as a clear example.
verify_b2_for_sun(3)