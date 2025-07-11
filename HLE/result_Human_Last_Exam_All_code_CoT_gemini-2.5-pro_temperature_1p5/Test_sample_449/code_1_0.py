import math

# Step 1: Define the starting point and calculate its distance from the origin.
x, y = 3000, 4000
r = math.sqrt(x**2 + y**2)

# Step 2: Define the constants needed for the calculation.
# Euler-Mascheroni constant
gamma = 0.57721566490153286

# Step 3: Calculate the terms in the formula.
# The formula for the probability P is: 1 - 1 / (ln(r) + gamma + 1.5 * ln(2))
ln_r = math.log(r)
ln_2 = math.log(2)

# The constant part of the denominator
constant_term = gamma + 1.5 * ln_2

# The full denominator
denominator = ln_r + constant_term

# The probability
prob = 1 - 1 / denominator

# Step 4: Print the numbers used in the calculation and the final result.
print("The probability P is given by the formula: P = 1 - a(1,0) / a(x0)")
print("We use the exact value a(1,0) = 2/pi")
print("We use the asymptotic formula a(x0) ≈ (2/pi) * (ln(|x0|) + γ + 1.5*ln(2))")
print("This simplifies to: P ≈ 1 - 1 / (ln(|x0|) + γ + 1.5*ln(2))\n")

print(f"Starting point x0 = ({x}, {y})")
print(f"Distance from origin |x0| = {r}")
print(f"ln(|x0|) = {ln_r:.4f}")
print(f"γ (Euler-Mascheroni constant) = {gamma:.4f}")
print(f"1.5 * ln(2) = {1.5 * ln_2:.4f}")
print(f"Denominator = ln(|x0|) + γ + 1.5*ln(2) = {ln_r:.4f} + {gamma:.4f} + {1.5*ln_2:.4f} = {denominator:.4f}")
print(f"\nThe probability is P ≈ 1 - 1 / {denominator:.4f}")
print(f"P ≈ {prob:.4f}")

# Final answer with two significant digits
approx_prob = f"{prob:.2g}"
print(f"\nThe approximate answer with two significant digits is: {approx_prob}")
