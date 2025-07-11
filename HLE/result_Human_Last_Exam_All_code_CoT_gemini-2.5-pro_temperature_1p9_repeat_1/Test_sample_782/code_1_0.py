# Plan: Compute 2A - 3B using the formula 2(A - B) - B.
# Let B_prime = -B. The calculation is equivalent to 2*(A + B_prime) + B_prime.
# The cost is measured in the number of field multiplications (M), where 1 Squaring (S) = 1M.

# Cost of primitive operations on twisted Edwards curves in extended coordinates
# based on the Explicit-Formulas Database (EFD).

# 1. Precomputation
# To work with extended coordinates, we often need the product of the affine coordinates.
# A is (x_A, y_A), B is (x_B, y_B). B_prime is (-x_B, y_B).
cost_T_A = 1  # To compute T_A = x_A * y_A
cost_T_B_prime = 1 # To compute T_B_prime = -x_B * y_B
cost_precomputation = cost_T_A + cost_T_B_prime
print(f"Step 1: Precomputation cost = {cost_precomputation}M")

# 2. First Addition: R1 = A + B_prime
# We add A and B_prime, both coming from affine coordinates.
# We can use their extended form where Z=1.
# The cost for adding two extended coordinate points with Z1=Z2=1 is 7M.
cost_first_addition = 7
print(f"Step 2: Cost of (A + B_prime) = {cost_first_addition}M")

# 3. Doubling: R2 = 2 * R1
# R1 is a general point in extended coordinates.
# The cost of doubling a general point is 8M.
cost_doubling = 8
print(f"Step 3: Cost of 2 * (A + B_prime) = {cost_doubling}M")

# 4. Second Addition: Result = R2 + B_prime
# This is a mixed addition of a general extended point (R2) and an affine point (B_prime).
# The general cost is 8M, but it's 7M if the product of affine coordinates (T_B_prime) is pre-computed.
# We already computed T_B_prime in Step 1.
cost_second_addition = 7
print(f"Step 4: Cost of adding B_prime = {cost_second_addition}M")

# Calculate the total cost
total_cost = cost_precomputation + cost_first_addition + cost_doubling + cost_second_addition

# Output the final equation
print("\nThe final cost is the sum of these steps:")
print(f"{cost_precomputation}M (Precomputation) + {cost_first_addition}M (First Addition) + {cost_doubling}M (Doubling) + {cost_second_addition}M (Second Addition) = {total_cost}M")
print(f"\nFinal Equation: {cost_precomputation} + {cost_first_addition} + {cost_doubling} + {cost_second_addition} = {total_cost}")
print(f"The smallest cost is {total_cost} multiplications.")
# Final answer format:
# <<<answer>>>
# where answer is the computed number.
# e.g. <<<24>>>
final_answer = total_cost
print(f'<<<{final_answer}>>>')
