import math

def calculate_fixed_point_coupling():
    """
    Calculates and prints the leading order expression for the fixed point
    coupling u* in phi^4 theory near d=4 dimensions.
    """
    # The beta function for the coupling u in d=4-epsilon dimensions is:
    # beta(u) = -epsilon * u + (3 / (16 * pi^2)) * u^2
    # The non-trivial fixed point u* is found by setting beta(u*) = 0.
    # -epsilon * u* + (3 / (16 * pi^2)) * (u*)^2 = 0
    # This gives u* = (16 * pi^2 / 3) * epsilon.

    # We will now calculate the numerical coefficient C = 16 * pi^2 / 3.
    
    numerator_coeff = 16
    denominator_coeff = 3
    
    # Calculate the coefficient
    coefficient = (numerator_coeff * math.pi**2) / denominator_coeff
    
    print("The leading order expression for the fixed point coupling u* is given by:")
    print("u* = C * epsilon")
    print("\nwhere the coefficient C is derived from the one-loop beta function.")
    print(f"The symbolic form of the coefficient is: C = ({numerator_coeff} * pi^2) / {denominator_coeff}")
    
    print("\nThe final expression with the calculated numerical coefficient is:")
    # The final print statement outputs each number in the final equation as requested.
    print(f"u* = ({numerator_coeff} * {math.pi}**2 / {denominator_coeff}) * epsilon = {coefficient} * epsilon")

calculate_fixed_point_coupling()
