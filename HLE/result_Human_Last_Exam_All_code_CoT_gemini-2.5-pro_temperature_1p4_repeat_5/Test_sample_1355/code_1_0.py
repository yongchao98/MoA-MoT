from fractions import Fraction

# Step 1: Define the values of the moments based on the derivation.
# The average conductance <g> is 1/2.
avg_g = Fraction(1, 2)

# The fourth moment of conductance <g^4> is 35/128.
g4 = Fraction(35, 128)

# Step 2: Calculate the ratio of the fourth moment to the average value.
ratio = g4 / avg_g

# Step 3: Print the results, showing each number in the final equation.
print(f"The average value of the dimensionless conductance is <g> = {avg_g}")
print(f"The fourth statistical moment of the dimensionless conductance is <g^4> = {g4}")
print("\nThe ratio between the fourth moment and the average value is <g^4> / <g>.")
print(f"The final equation is: {g4} / {avg_g} = {ratio}")
print(f"The numerical value of the ratio is approximately {float(ratio):.5f}")
