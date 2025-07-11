import math

# --- Step 1: Define constants and parameters ---
# The starting position is (3000, 4000).
# The distance (r) from the origin is sqrt(3000^2 + 4000^2) = 5000.
r = 5000.0

# The Euler-Mascheroni constant (gamma_E)
gamma_E = 0.5772156649

# --- Step 2: Calculate the terms in the formula ---
# The formula for the probability is log(r) / (log(r) + gamma_E + log(8))
log_r = math.log(r)
log_8 = math.log(8)

# --- Step 3: Compute the final probability ---
numerator = log_r
denominator = log_r + gamma_E + log_8
probability = numerator / denominator

# --- Step 4: Print the results ---
print("This script calculates the approximate probability that a 2D random walk,")
print("conditioned to avoid the origin and starting at (3000, 4000),")
print("will never enter the set of the four neighbors of the origin.")
print("\n--- Calculation ---")
print(f"The distance from the origin is r = {r}")
print(f"The formula for the probability is: log(r) / (log(r) + gamma_E + log(8))")
print("\nPlugging in the values:")
# We print each number in the final equation.
print(f"P ≈ {numerator:.4f} / ({log_r:.4f} + {gamma_E:.4f} + {log_8:.4f})")
print(f"P ≈ {numerator:.4f} / {denominator:.4f}")
print(f"P ≈ {probability:.4f}")

# --- Step 5: Final Answer ---
# The question asks for the approximate answer with two significant digits.
final_answer = round(probability, 2)
print(f"\nThe approximate answer with two significant digits is: {final_answer}")
<<<0.76>>>