# The genus of the Riemann surface C.
g = 3

# The symmetric power d.
d = 15

print("The problem is to find the smallest and largest possible ranks of the Neron-Severi group of X = C^(15),")
print(f"where C is a Riemann surface of genus g = {g}.\n")

print("The rank of the Neron-Severi group of the d-th symmetric power of C, rho(C^(d)), is given by the formula:")
print("rho(C^(d)) = rho(J(C)) + 1")
print("where J(C) is the Jacobian of C. This formula holds for d >= g-1.\n")
print(f"In this case, g = {g} and d = {d}, so {d} >= {g-1}, and the formula applies.\n")

print(f"The rank of the Neron-Severi group of the Jacobian, rho(J(C)), for a curve of genus g={g} is bounded:")
print(f"1 <= rho(J(C)) <= g^2, which means 1 <= rho(J(C)) <= {g**2}.\n")

# --- Smallest Rank ---
min_rho_JC = 1
smallest_rank = min_rho_JC + 1
print("1. Smallest Rank:")
print(f"This occurs for a generic curve C, where rho(J(C)) is minimal, i.e., rho(J(C)) = {min_rho_JC}.")
print("The smallest possible rank of NS(X) is therefore:")
print(f"rho(J(C))_min + 1 = {min_rho_JC} + 1 = {smallest_rank}\n")

# --- Largest Rank ---
max_rho_JC = g**2
largest_rank = max_rho_JC + 1
print("2. Largest Rank:")
print(f"This occurs for a special curve C (like the Fermat quartic), where rho(J(C)) is maximal, i.e., rho(J(C)) = g^2 = {max_rho_JC}.")
print("The largest possible rank of NS(X) is therefore:")
print(f"rho(J(C))_max + 1 = {max_rho_JC} + 1 = {largest_rank}\n")
