import sys

def verify_betti_number_formula():
    """
    This function investigates the claim that for G = SU(n), the second
    Betti number b_2 of any coadjoint orbit is given by n - 1.

    We will test this claim with a counterexample.
    """
    
    # Consider the group G = SU(3), so n = 3.
    n = 3
    print(f"We are examining the case G = SU(n) where n = {n}.")

    # The formula in question is b_2 = n - 1.
    # Let's calculate the predicted value.
    predicted_b2 = n - 1
    
    print(f"The formula b_2 = n - 1 predicts the second Betti number to be: {n} - 1 = {predicted_b2}")

    # Now, we consider a specific coadjoint orbit for SU(3). Not all orbits are
    # the full flag manifold SU(3)/T. If we choose a 'singular' lambda,
    # we can get smaller orbits. One such orbit is diffeomorphic to the
    # complex projective plane, CP^2.
    orbit_name = "the complex projective plane (CP^2)"
    
    # The Betti numbers of CP^2 are well-known in topology.
    # b_0 = 1, b_1 = 0, b_2 = 1, b_3 = 0, b_4 = 1.
    actual_b2 = 1
    
    print(f"Let's consider a singular coadjoint orbit, which is {orbit_name}.")
    print(f"The actual second Betti number of {orbit_name} is {actual_b2}.")

    # Finally, we compare the prediction with the actual value.
    print("\nConclusion for part (b):")
    if predicted_b2 == actual_b2:
        print(f"The values are equal ({predicted_b2} == {actual_b2}). The formula holds in this case.")
    else:
        print(f"The values are not equal ({predicted_b2} != {actual_b2}).")
        print("Therefore, the statement that b_2 is *always* n - 1 is false.")


# Execute the verification
if __name__ == '__main__':
    verify_betti_number_formula()