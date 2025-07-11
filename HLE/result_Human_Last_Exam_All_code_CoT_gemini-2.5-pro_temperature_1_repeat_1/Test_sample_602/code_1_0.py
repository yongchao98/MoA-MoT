# All calculations are done symbolically, leading to a final formula for l(n).
# This Python code serves to present the final derived formula and its components as requested.

# The derivation shows that l(n) can be expressed as a function of n.
# The formula is: l(n) = (c1*n^2 + c2 - (c3*n - c4) * sqrt(c5*n^2 - c6*n + c7)) / n^2

# Define the coefficients of the derived formula.
coefficients = {
    'c1': 2,  # Coefficient for n^2 in the numerator
    'c2': 2,  # Constant term in the numerator
    'c3': 4,  # Coefficient for n inside the parenthesis before the square root
    'c4': 2,  # Constant term inside the parenthesis
    'c5': 1,  # Coefficient for n^2 inside the square root
    'c6': 1,  # Coefficient for n inside the square root
    'c7': 1,  # Constant term inside the square root
}

# Print the final equation for l(n)
print("The exact value of l(n) is given by the equation:")
print(f"l(n) = ({coefficients['c1']}*n^2 + {coefficients['c2']} - ({coefficients['c3']}*n - {coefficients['c4']}) * sqrt({coefficients['c5']}*n^2 - {coefficients['c6']}*n + {coefficients['c7']})) / n^2")

# Print each numerical coefficient as requested.
print("\nThe numbers in the final equation are:")
for name, value in coefficients.items():
    print(f"  {name}: {value}")
