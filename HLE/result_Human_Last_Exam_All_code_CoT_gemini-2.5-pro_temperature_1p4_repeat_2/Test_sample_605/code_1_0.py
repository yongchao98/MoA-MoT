def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a Calabi-Yau link.
    
    The invariant is defined as d - sum(w_i), where d is the weighted
    degree of the polynomial and w_i are the weights of the ambient space.
    """
    
    # Weights for the variables (z1, z2, z3, z4, z5)
    weights = [22, 29, 49, 50, 75]
    
    # The weighted degree 'd' of the quasi-homogeneous polynomial.
    # As determined by analyzing the monomials, the consistent degree is 225.
    degree = 225
    
    # Calculate the sum of the weights
    sum_of_weights = sum(weights)
    
    # Calculate the invariant
    invariant = degree - sum_of_weights
    
    # Construct the weight sum string for the final equation
    weights_str = " + ".join(map(str, weights))
    
    print("The Crawley-Nordström invariant is calculated as the degree (d) minus the sum of the weights (Σwᵢ).")
    print(f"For the given link:")
    print(f"Degree (d) = {degree}")
    print(f"Weights (wᵢ) = {weights}")
    print(f"The calculation is:")
    print(f"{degree} - ({weights_str}) = {invariant}")
    
    # Final answer in the specified format
    print(f"\n<<<The final answer is {invariant}>>>")

calculate_crawley_nordstrom_invariant()