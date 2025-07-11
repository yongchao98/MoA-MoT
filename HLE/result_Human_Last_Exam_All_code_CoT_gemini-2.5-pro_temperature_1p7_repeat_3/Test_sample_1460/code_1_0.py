def solve_knot_problem():
    """
    Solves the knot theory problem by identifying the components of the link
    and determining the knot type of the specified component.
    """

    # The braid is given by beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^{-1} in the braid group B_5.
    braid_word = "s1^2 * s2^2 * s3 * s4^-1"
    
    # Step 1: Find the permutation of the braid to identify the link components.
    # pi(s1^2) is the identity permutation.
    # pi(s2^2) is the identity permutation.
    # pi(s3) is the transposition (3, 4).
    # pi(s4^-1) is the transposition (4, 5).
    # The total permutation is (3, 4) * (4, 5) = (3, 4, 5).
    # The cycle decomposition is (1)(2)(3, 4, 5).
    # This means the link has 3 components: C1={1}, C2={2}, C3={3, 4, 5}.
    
    # Step 2: The problem states two components are unknots. These are C1 and C2.
    # We need to find the knot type of the third component, C3.
    
    # Step 3: The knot type of C3 is determined by the braiding among its own strands {3, 4, 5}.
    # The relevant part of the braid word is sigma_3 * sigma_4^{-1}.
    
    # Step 4: We analyze this as a 3-strand braid. Relabeling strands {3,4,5} to {1,2,3},
    # sigma_3 becomes sigma_1 and sigma_4 becomes sigma_2.
    # The equivalent braid in B_3 is sigma_1 * sigma_2^{-1}.
    
    # Step 5: The closure of the braid sigma_1 * sigma_2^{-1} is a well-known knot.
    # This knot is the figure-8 knot, also known as 4_1.
    
    knot_type = "Figure-8"
    
    print(f"The braid is beta = {braid_word} in B_5.")
    print("The closure of the braid has 3 connected components.")
    print("Given that two components are unknots, the third component is formed by strands {3, 4, 5}.")
    print("The braiding on these strands corresponds to the 3-braid sigma_1 * sigma_2^{-1}.")
    print(f"The closure of this braid is the {knot_type} knot.")

solve_knot_problem()