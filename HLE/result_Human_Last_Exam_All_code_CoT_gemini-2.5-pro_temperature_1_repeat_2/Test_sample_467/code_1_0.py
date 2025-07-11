import sympy

# Define the variable z
z = sympy.symbols('z')

# Define the numerator and denominator of the Gauss map g(z) = z / (z^3 + 2)
P = z
Q = z**3 + 2

# 1. Calculate the degree 'd' of the Gauss map
# The degree of a rational function is the maximum of the degrees of the numerator and denominator.
deg_P = sympy.degree(P, gen=z)
deg_Q = sympy.degree(Q, gen=z)
d = max(deg_P, deg_Q)

# 2. Calculate the number of ends 'k'
# The number of ends is the number of poles of the Gauss map, which is the degree of the denominator.
# We confirm the roots are distinct by checking that the gcd of Q and its derivative is a constant.
# Q_prime = sympy.diff(Q, z) -> 3*z**2. The only root is z=0. Q(0) = 2 != 0. So roots are distinct.
k = sympy.degree(Q, gen=z)

# 3. Apply the López-Pérez formula for the Morse index: Index = 2d - k - 1
index = 2 * d - k - 1

# Print the results step-by-step
print(f"The Gauss map is g(z) = z / (z^3 + 2).")
print(f"The degree of the Gauss map is d = max(degree(z), degree(z^3 + 2)) = {d}.")
print(f"The number of ends 'k' is the number of poles, which is the degree of the denominator = {k}.")
print("\nUsing the formula: Index = 2*d - k - 1")
print(f"Index = 2 * {d} - {k} - 1 = {index}")

# Final Answer
print(f"\nThe Morse index of M is {index}.")