import numpy as np

# The infinite product prod_{n=0 to inf}(1 - e^(-(2n+1)pi)) has been shown to
# have the closed-form expression (2 / e^pi)^(1/24).
# This script calculates the numerical value of this expression.

# The final equation for the product P is: P = (A / e^B)^C
A = 2
B = np.pi
C = 1/24

print("The closed-form expression for the infinite product is (2 / e^pi)^(1/24).")
print("This can be written in the form: P = (A / e^B)^C")
print("\nThe numbers in this final equation are:")
print(f"A = {A}")
print(f"B = pi = {B}")
print(f"C = 1/24 = {C}")

# Calculate the final value
result = (A / np.exp(B))**C

print(f"\nThe numerical value of the expression is:")
print(f"P = ({A} / e^{B})^{C} = {result}")
