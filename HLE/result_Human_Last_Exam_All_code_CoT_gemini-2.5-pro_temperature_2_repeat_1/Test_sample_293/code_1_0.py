import numpy as np

def solve_integral():
    """
    Solves the integral by simplifying the integrand and using a known mathematical result.
    """
    
    print("The integral is I = integral from 0 to inf of sum_{n=1 to inf} log(cos(x/2^n)) dx.")
    print("-" * 20)

    # Step 1: Simplify the sum in the integrand.
    # The sum can be written as a logarithm of an infinite product:
    # sum_{n=1 to inf} log(cos(x/2^n)) = log( product_{n=1 to inf} cos(x/2^n) )
    
    # Step 2: Use the known product identity: product_{n=1 to inf} cos(x/2^n) = sin(x) / x.
    # Let's demonstrate this identity numerically for x=1.5.
    x_test = 1.5
    
    # Calculate the product side
    product_val = 1.0
    for n in range(1, 50):  # Using 50 terms for good approximation
        product_val *= np.cos(x_test / (2**n))
        
    # Calculate the sin(x)/x side
    sinc_val = np.sin(x_test) / x_test
    
    print(f"Step 2: Numerical check of the identity at x = {x_test}")
    print(f"Infinite product approximation: {product_val}")
    print(f"sin(x)/x value:               {sinc_val}")
    print("The values are very close, which supports the identity.")
    print("-" * 20)
    
    # Step 3: The integral becomes I = integral from 0 to inf of log(sin(x)/x) dx.
    print("Step 3: The integral simplifies to the Lobachevsky integral.")
    print("I = integral from 0 to inf of log(sin(x)/x) dx")
    print("-" * 20)

    # Step 4: This is a known integral whose value is -pi/2.
    print("Step 4: The value of this standard integral is -pi / 2.")
    
    pi = np.pi
    numerator = -pi
    denominator = 2
    
    result = numerator / denominator
    
    print("The final calculation is:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")
    print(f"Result: {numerator} / {denominator} = {result}")

solve_integral()