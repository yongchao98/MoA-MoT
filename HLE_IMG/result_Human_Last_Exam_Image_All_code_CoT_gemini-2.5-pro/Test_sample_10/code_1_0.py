# Data points for the inversion barrier (B) vs. molecule index (n)
# n=1 corresponds to the first molecule, n=2 to the second, and so on.
b1 = 10  # Barrier for n=1
b2 = 49  # Barrier for n=2

# Index of the molecule for which we want to predict the barrier
n3 = 3

# Explanation of the model
print("This script predicts the inversion barrier for the third molecule in the series.")
print("The prediction is based on a quadratic model derived from the first two molecules.")
print("The model assumes the barrier B is related to the molecule index n by the formula: B(n) = a * n^2 + c.")
print("-" * 50)

# We solve the system of equations for the coefficients 'a' and 'c':
# For n=1: a * 1^2 + c = 10
# For n=2: a * 2^2 + c = 49
#
# Solving (4a + c = 49) - (a + c = 10) gives 3a = 39, so a = 13.
# Substituting 'a' back gives 13 + c = 10, so c = -3.
a = (b2 - b1) / (2**2 - 1**2)
c = b1 - a * 1**2

# Ensure coefficients are integers for clean output
a = int(a)
c = int(c)

# Calculate the predicted barrier for n=3
predicted_barrier = a * (n3**2) + c

# Print the derivation and final result, showing each number in the equation.
print("The derived model is: B(n) = {} * n^2 + ({})".format(a, c))
print(f"\nCalculating the barrier for the third molecule (n={n3}):")
print(f"Barrier = {a} * {n3}^2 + ({c})")
print(f"        = {a} * {n3**2} + ({c})")
print(f"        = {a * n3**2} + ({c})")
print(f"        = {predicted_barrier}")
print(f"\nThe predicted inversion barrier is {predicted_barrier} kcal/mol.")
<<<114>>>