import numpy as np

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau link.
    """
    # The weights for the variables (z_1, z_2, z_3, z_4, z_5)
    weights = np.array([22, 29, 49, 50, 75])
    
    # Exponents for each term in the polynomial W.
    # W = z_1^8z_3 + z_1^4z_2^3z_3 + z_1z_2^7 + z_1z_2z_3z_4z_5 + z_2z_3^4 + z_4^3z_5 + z_5^3
    # Note: As discussed, there is likely a typo in the second term.
    # We will assume the degree d is determined by the Calabi-Yau condition (d = sum of weights).
    
    # 1. Calculate the sum of the weights
    sum_of_weights = np.sum(weights)
    
    # 2. Determine the degree d of the polynomial
    # For a Calabi-Yau hypersurface, the degree 'd' equals the sum of the weights.
    d = sum_of_weights
    
    # 3. Calculate the Crawley-Nordström invariant
    # Formula: inv(W) = sum(w_i) - 2 * d
    invariant = sum_of_weights - 2 * d
    
    # Print the explanation and the final equation
    print("The Crawley-Nordström invariant is calculated using the formula:")
    print("Invariant = (sum of weights) - 2 * (degree d)")
    print("\nGiven values:")
    print(f"Weights = {list(weights)}")
    print(f"Sum of weights = {sum_of_weights}")
    print(f"Degree (d) = {d}")
    print("\nCalculation:")
    # The user requested to see each number in the final equation.
    print(f"{sum_of_weights} - 2 * {d} = {invariant}")

calculate_crawley_nordstrom_invariant()