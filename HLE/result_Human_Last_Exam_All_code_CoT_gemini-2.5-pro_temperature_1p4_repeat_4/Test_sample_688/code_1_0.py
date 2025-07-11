import math
from fractions import Fraction

# The user can change this value for the number of nodes.
n = 4

if not isinstance(n, int) or n < 2:
    print(f"Error: The number of nodes n must be an integer greater than or equal to 2. Got n = {n}")
else:
    # The prefactor c_n is given by the formula -(n-1)/n!
    
    # Calculate numerator and denominator of the unsimplified fraction
    numerator_val = -(n - 1)
    denominator_val = math.factorial(n)
    
    # Use the Fraction class for an exact fractional representation
    result_fraction = Fraction(numerator_val, denominator_val)
    
    # Print the explanation and the calculation steps
    print(f"The system-independent prefactor c_n for the fully f-connected diagram with n nodes is given by the formula:")
    print("c_n = -(n - 1) / n!")
    print(f"\nFor the specific case of n = {n}:")
    
    # Output the equation with numbers plugged in
    print(f"c_{n} = -({n} - 1) / {n}!")
    
    # Output the intermediate calculation before simplification
    print(f"c_{n} = {numerator_val} / {denominator_val}")
    
    # Output the final simplified fraction and its decimal value
    print(f"The simplified prefactor is: c_{n} = {result_fraction}")
    print(f"As a decimal: c_{n} = {float(result_fraction)}")
