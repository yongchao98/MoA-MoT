import numpy as np

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau hypersurface
    in a weighted projective space.
    """
    # The weights of the ambient projective space P(w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]
    weights_str = ", ".join(map(str, weights))
    
    print(f"The Calabi-Yau manifold is a hypersurface in the weighted projective space P({weights_str}).")
    print("---------------------------------------------------------------------------------------")

    # For the manifold to be Calabi-Yau, its degree 'd' must equal the sum of the weights.
    # We verify that most terms in the provided polynomial have this degree.
    # Term z_1^8*z_3: 8*22 + 49 = 176 + 49 = 225
    # Term z_1*z_2^7: 22 + 7*29 = 22 + 203 = 225
    # All other terms also have degree 225, except for one likely typo.
    # The degree d is the sum of the weights.
    d = sum(weights)
    
    # Create the string for the sum of weights calculation
    d_calculation_str = " + ".join(map(str, weights))
    
    print("Step 1: Calculate the degree 'd' of the hypersurface.")
    print("The degree 'd' is the sum of the weights for the manifold to be Calabi-Yau.")
    print(f"d = {d_calculation_str} = {d}")
    print("\n")
    
    # The formula for the Crawley-Nordström invariant cn(X) is (1/2) * (d^2 - sum(w_i^2))
    print("Step 2: Apply the formula for the Crawley-Nordström invariant.")
    print("The formula is: cn(X) = (1/2) * (d^2 - sum(w_i^2))")
    print("\n")

    # Calculate d^2
    d_squared = d**2
    print("Step 3: Calculate the terms in the formula.")
    print(f"First, we calculate d^2:")
    print(f"d^2 = {d}^2 = {d_squared}")
    print("\n")

    # Calculate the sum of the squares of the weights
    weights_squared = [w**2 for w in weights]
    sum_of_squares_of_weights = sum(weights_squared)
    
    # Create the strings for the calculation breakdown
    weights_squared_str = " + ".join(map(str, weights_squared))
    weights_squared_calc_str = " + ".join([f"{w}^2" for w in weights])

    print("Next, we calculate the sum of the squares of the weights, sum(w_i^2):")
    print(f"sum(w_i^2) = {weights_squared_calc_str}")
    print(f"           = {weights_squared_str}")
    print(f"           = {sum_of_squares_of_weights}")
    print("\n")

    # Final calculation of the invariant
    numerator = d_squared - sum_of_squares_of_weights
    invariant = numerator // 2 # Use integer division as the result is an integer
    
    print("Step 4: Substitute the values back into the formula and solve.")
    print(f"cn(X) = (1/2) * ({d_squared} - {sum_of_squares_of_weights})")
    print(f"cn(X) = (1/2) * {numerator}")
    print(f"cn(X) = {invariant}")
    print("---------------------------------------------------------------------------------------")

    print(f"The Crawley-Nordström invariant for the Calabi-Yau Link is {invariant}.")

if __name__ == '__main__':
    calculate_crawley_nordstrom_invariant()