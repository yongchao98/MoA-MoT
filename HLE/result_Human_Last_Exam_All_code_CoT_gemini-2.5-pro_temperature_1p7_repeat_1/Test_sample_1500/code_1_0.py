def check_betti_number_claim(n):
    """
    This function checks the claim from part (b) of the question.
    The claim is that for G = SU(n), the second Betti number b_2 of any
    coadjoint orbit is n - 1.

    We demonstrate this is false with a counterexample. There exist
    coadjoint orbits of SU(n) that are diffeomorphic to the complex
    projective space CP^{n-1}.

    The Betti numbers of CP^{k} are b_{2i} = 1 for 0 <= i <= k and 0 otherwise.
    So, the second Betti number b_2(CP^{n-1}) is 1.

    This function compares this value to the claimed value n - 1.
    """
    if not isinstance(n, int) or n < 2:
        print("Please provide an integer n >= 2 for SU(n).")
        return

    # For the coadjoint orbit diffeomorphic to CP^{n-1}
    b2_of_orbit = 1
    
    # The value claimed in the question
    claimed_b2 = n - 1
    
    print(f"For G = SU({n}):")
    print(f"Consider the coadjoint orbit O_lambda which is diffeomorphic to CP^{n-1}.")
    print(f"The second Betti number of this orbit is b_2(CP^{n-1}) = {b2_of_orbit}.")
    print(f"The question claims b_2 is always {n} - 1 = {claimed_b2}.")
    
    if b2_of_orbit == claimed_b2:
        print(f"For n = {n}, the values match: {b2_of_orbit} == {claimed_b2}.")
        print("This holds for n=2.")
    else:
        print(f"For n = {n}, the values do NOT match: {b2_of_orbit} != {claimed_b2}.")
        print("This provides a counterexample to the claim for n > 2.")

# Run the check for a specific case, e.g., n = 3, as discussed in the explanation.
check_betti_number_claim(3)

print("\n---")
# Run the check for n=2, where the claim happens to be true.
check_betti_number_claim(2)