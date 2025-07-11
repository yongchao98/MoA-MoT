import math

# Assume n is a given positive integer.
# For demonstration, we'll set n to a value.
# The user can change this value to their desired n.
n = 10

# The formula for the maximum number of eigenvalues > 2 is n - floor(3 / (n + 1)).
# The following code calculates this value and prints the components of the
# final equation as per the instructions.

if not isinstance(n, int) or n <= 0:
    print("Error: n must be a positive integer.")
else:
    # These are the constant numbers in the derived formula.
    numerator = 3
    denominator_add = 1
    
    # In Python, integer division `//` on positive numbers is equivalent to floor().
    result = n - (numerator // (n + denominator_add))
    
    print(f"For n = {n}:")
    print(f"The maximum number of eigenvalues greater than 2 is {result}.")
    print("\nThis result is calculated from the expression: n - floor(3 / (n + 1))")
    print("The final equation with its numerical components is:")
    
    # This line prints out each number in the final equation as requested.
    print(f"{n} - floor({numerator} / ({n} + {denominator_add}))")
