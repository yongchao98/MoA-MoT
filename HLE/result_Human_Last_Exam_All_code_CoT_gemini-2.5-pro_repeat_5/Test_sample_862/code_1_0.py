import numpy as np

# Define theta and k from our analytical solution
theta = 2 * np.pi / 3
k = 0.5

# The transcendental equation derived for the odd eigenfunction:
# sin(k*(2*theta - pi)) = 0.5 * sin(k*pi)
lhs = np.sin(k * (2 * theta - np.pi))
rhs = 0.5 * np.sin(k * np.pi)

print(f"For theta = {theta:.4f} and k = {k}")
print(f"LHS of the equation: sin(k*(2*theta - pi)) = {lhs:.6f}")
print(f"RHS of the equation: 0.5 * sin(k*pi) = {rhs:.6f}")

if np.isclose(lhs, rhs):
    print("\nThe equation is satisfied.")
    C = 1 / k**2
    print(f"The eigenvalue is k^2 = {k**2}, which gives the constant C = 1/k^2 = {C}")
else:
    print("\nThe equation is NOT satisfied.")

# Print the final equation values for the user
print("\nFinal calculation:")
print(f"C = 1 / ({k})**2 = {1/k**2}")
