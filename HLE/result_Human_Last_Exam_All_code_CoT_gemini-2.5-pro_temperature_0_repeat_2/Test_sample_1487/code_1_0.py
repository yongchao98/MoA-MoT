import math

# Step 1: Define the known infinite sum from the Basel problem.
# S_1 = sum_{k=1 to inf} (1/k^2) = pi^2 / 6
basel_sum = (math.pi**2) / 6
print(f"The Basel sum is S₁ = Σ(1/k²) from k=1 to ∞ = π²/6 ≈ {basel_sum}")

# Step 2: Calculate the sum needed for our problem.
# S_2 = sum_{i=1 to inf} (1/(i+1)^2) = sum_{k=2 to inf} (1/k^2) = S_1 - 1
sum_from_2 = basel_sum - 1
print(f"The required sum is S₂ = Σ(1/(i+1)²) from i=1 to ∞ = S₁ - 1 ≈ {sum_from_2}")

# Step 3: Calculate the squared norm of the vector α.
# ||α||² = (1/2) * S_2
norm_alpha_sq = 0.5 * sum_from_2
print(f"The squared norm is ||α||² = (1/2) * S₂ ≈ {norm_alpha_sq}")

# Step 4: Evaluate the final expression.
# Expression = (2 * ||α||²) / (π²/6 - 1) + 10^15
# Note that π²/6 - 1 is just S_2.
# Expression = (2 * ||α||²) / S₂ + 10^15
numerator = 2 * norm_alpha_sq
denominator = sum_from_2
constant = 10**15

# The calculation simplifies because numerator = 2 * (0.5 * S_2) = S_2.
# So, (numerator / denominator) = 1.
final_value = (numerator / denominator) + constant

print("\nCalculating the final expression:")
print(f"Expression = (2 * ||α||²) / (π²/6 - 1) + 10¹⁵")
print(f"             = (2 * {norm_alpha_sq}) / ({denominator}) + {constant}")
print(f"             = ({numerator}) / ({denominator}) + {constant}")
print(f"             = {numerator / denominator} + {constant}")
print(f"Final Value  = {final_value}")

# The final result is an integer, so we can print it as such.
print(f"\nThe final result is: {int(final_value)}")