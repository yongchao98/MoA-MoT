# This script demonstrates the counterexample for question (b).
# It checks if the second Betti number b_2 of a coadjoint orbit of SU(n)
# is always given by the formula n - 1.
# We will test this for n=4 using known results from topology.

def check_betti_number_formula(n):
    """
    Checks the formula b_2 = n - 1 for two different SU(n) coadjoint orbits.
    """
    predicted_b2 = n - 1

    print(f"Analysis for G = SU({n}):")
    print(f"The question is whether the second Betti number b_2 is always {n} - 1 = {predicted_b2}.")
    print("-" * 50)

    # Case 1: Regular Orbit (Full Flag Manifold)
    # The orbit is SU(n)/T, where T is the maximal torus.
    # The second Betti number b_2(SU(n)/T) equals the rank of SU(n), which is n-1.
    case1_name = f"Regular Orbit (Full Flag Manifold SU({n})/T)"
    case1_actual_b2 = n - 1
    print(f"Case 1: {case1_name}")
    print(f"The actual second Betti number is rank(SU({n})) = {case1_actual_b2}.")
    print(f"The formula b_2 = {n} - 1 gives {predicted_b2}.")
    if case1_actual_b2 == predicted_b2:
        print("Result: The formula holds for this case.")
    else:
        print("Result: The formula fails for this case.")
    print("-" * 50)

    # Case 2: Singular Orbit (Grassmannian Gr(2, n))
    # This orbit is SU(n)/S(U(2) x U(n-2)).
    # It is a known result that b_2(Gr(k, n)) = 1 for any k. We take k=2.
    if n > 2:
        k = 2
        case2_name = f"Singular Orbit (Grassmannian Gr({k}, {n}))"
        case2_actual_b2 = 1
        print(f"Case 2: {case2_name}")
        print(f"The actual second Betti number is {case2_actual_b2}.")
        print(f"The formula b_2 = {n} - 1 gives {predicted_b2}.")
        if case2_actual_b2 == predicted_b2:
            print("Result: The formula holds for this case.")
        else:
            print("Result: The formula fails for this case.")
        print("-" * 50)

# We run the check for n=4, which provides a clear counterexample.
check_betti_number_formula(4)
print("\nConclusion: Since we found a case where the formula fails, the statement is false.")
