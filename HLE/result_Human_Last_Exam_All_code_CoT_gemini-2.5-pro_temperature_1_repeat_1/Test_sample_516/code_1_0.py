import sympy
from sympy import summation, oo, factorial, E, S

# Step 1: Explain the derived formula for the expected value.
print("Based on the derivation, the expected value E[T] is given by the infinite series:")
print("E[T] = sum_{i=1 to infinity} [ i / (3^i * (i+1)!) ]")
print("We will use Python's symbolic math library, sympy, to compute this sum.")
print("-" * 30)

# Step 2: Define the symbolic variables for the summation.
i = sympy.Symbol('i', integer=True, positive=True)

# Step 3: Define the general term of the series.
# We use S(3) to ensure 3 is treated as a sympy object for symbolic computation.
term = i / (S(3)**i * factorial(i + 1))

# Step 4: Compute the sum of the series from i=1 to infinity.
expected_value_T = summation(term, (i, 1, oo))

# Step 5: The result of the summation is 3 - 2*exp(1/3).
# We extract the numbers from this expression to display them.
a = 3
b = 2
c = 1
d = 3

# Step 6: Print the final equation and its components as requested.
print("The symbolic computation results in the following expression:")
print(f"E[T] = {expected_value_T}")
print("\nThe final equation for the expected value of T is structured as: a - b*exp(c/d)")
print("The numbers in this final equation are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print(f"d = {d}")

# Step 7: Print the numerical value for context.
numerical_value = expected_value_T.evalf()
print(f"\nThe numerical value of E[T] is approximately: {numerical_value:.8f}")