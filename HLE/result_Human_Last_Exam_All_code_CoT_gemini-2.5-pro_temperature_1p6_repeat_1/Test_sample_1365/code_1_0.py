import math

# Define the parameters of the problem
# n: the total number of experts
n = 11
# c: the number of mistakes an expert can make before being removed
c = 4

# --- Calculation based on the derived formula ---
# Bound for M1 (mistakes when the true expert is wrong)
m1_bound = c - 1

# Bound for M2 (mistakes when the true expert is right)
m2_bound = c * (n - 1) / 2

# Total bound
total_bound = m1_bound + m2_bound

# --- Output the derivation with the given numbers ---
print("Deriving the upper bound for the algorithm's mistakes (M):")
print(f"Given n = {n} experts and a mistake limit c = {c}.")
print("\nThe total mistakes M <= M1 (true expert is wrong) + M2 (true expert is right).")

print("\n1. Bounding M1:")
print(f"   The true expert makes at most c - 1 mistakes.")
print(f"   M1 <= {c} - 1 = {int(m1_bound)}")

print("\n2. Bounding M2:")
print(f"   Each M2-type mistake requires at least 2 non-true experts to be wrong.")
print(f"   The total mistake budget for the n - 1 non-true experts is c * (n - 1).")
print(f"   Mistake budget = {c} * ({n} - 1) = {c * (n - 1)}")
print(f"   Therefore, 2 * M2 <= c * (n - 1).")
print(f"   M2 <= ({c} * ({n} - 1)) / 2 = {int(m2_bound)}")

print("\n3. Final Equation and Total Bound:")
print(f"   M <= (c - 1) + (c * (n - 1)) / 2")
print(f"   Substituting the values: M <= ({c} - 1) + ({c} * ({n} - 1)) / 2")
print(f"   M <= {int(m1_bound)} + {int(m2_bound)}")
print(f"   M <= {int(total_bound)}")