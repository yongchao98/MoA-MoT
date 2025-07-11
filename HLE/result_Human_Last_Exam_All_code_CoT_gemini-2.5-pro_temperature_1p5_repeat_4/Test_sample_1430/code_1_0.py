# The genus of the Riemann surface C
g = 3

# The rank of the Neron-Severi group of X = C^(15) is given by
# rho(X) = 1 + rho(J(C)), where J(C) is the Jacobian of C.
# We need to find the minimum and maximum possible values for rho(J(C)).

# For a g-dimensional Jacobian J(C), the rank of its Neron-Severi group
# is bounded by 1 <= rho(J(C)) <= g^2.

# Smallest possible rank of NS(J(C)) for a generic curve
rho_J_min = 1

# Largest possible rank of NS(J(C)) for a special curve with maximal Complex Multiplication
rho_J_max = g**2

# Calculate the smallest possible rank for NS(X)
rho_X_min = 1 + rho_J_min

# Calculate the largest possible rank for NS(X)
rho_X_max = 1 + rho_J_max

# Print the results, showing the full equations
print(f"The genus of the curve C is g = {g}.")
print("The rank of the Neron-Severi group of X=C^(15) is rho(X) = 1 + rho(J(C)).")
print("\n--- Smallest Possible Rank ---")
print("For a generic curve, the rank of the Jacobian's Neron-Severi group is 1.")
print(f"Smallest rank of NS(X) = 1 + rho(J(C))_min = 1 + {rho_J_min} = {rho_X_min}")

print("\n--- Largest Possible Rank ---")
print("For a special curve, the rank of the Jacobian's Neron-Severi group is g^2.")
print(f"Largest rank of NS(X) = 1 + rho(J(C))_max = 1 + {g}^2 = 1 + {rho_J_max} = {rho_X_max}")