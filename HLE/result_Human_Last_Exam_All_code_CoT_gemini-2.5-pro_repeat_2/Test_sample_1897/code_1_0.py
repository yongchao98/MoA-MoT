import math

def solve_infinite_product():
    """
    Calculates and prints the closed-form expression for the infinite product
    prod_{n=0 to inf} (1 - exp(-(2n+1)pi)).
    """
    
    # The closed form of the product is a known result from number theory.
    # It is equal to 2 raised to the power of -1/8.
    
    # Define the components of the final equation string
    product_str = "prod_{n=0 to inf}(1 - e^(-(2*n+1)*pi))"
    expression_str = "2**(-1/8)"
    
    # Calculate the numerical value
    value = 2**(-1/8)
    
    # Print the equation
    print(f"The closed form of the infinite product is:")
    print(f"{product_str} = {expression_str}")
    
    # Print the numerical value
    print(f"\nNumerical value:")
    print(value)

solve_infinite_product()