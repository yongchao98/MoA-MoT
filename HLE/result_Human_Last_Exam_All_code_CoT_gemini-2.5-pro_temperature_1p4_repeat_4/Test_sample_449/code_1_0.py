import math

# Step 1: Define constants
pi = math.pi
gamma = math.gamma(1) # Euler-Mascheroni constant
x0_norm = 5000

# Step 2: Calculate the terms of the final expression
# P = 1 - (4 - pi) / (2*log(5000) + 2*gamma + log(8) + 4)
numerator = 4 - pi
log_x0 = math.log(x0_norm)
log_8 = math.log(8)

denominator = 2 * log_x0 + 2 * gamma + log_8 + 4

# Step 3: Calculate the probability
prob_hit = numerator / denominator
prob_never_hit = 1 - prob_hit

# Step 4: Print the equation and the result
print("The probability P is given by the formula:")
print(f"P = 1 - (4 - pi) / (2*log({x0_norm}) + 2*gamma + log(8) + 4)")
print("\nSubstituting the values:")
print(f"4 - pi = {numerator:.4f}")
print(f"2*log({x0_norm}) = {2*log_x0:.4f}")
print(f"2*gamma = {2*gamma:.4f}")
print(f"log(8) = {log_8:.4f}")
print(f"Denominator = {2*log_x0:.4f} + {2*gamma:.4f} + {log_8:.4f} + 4 = {denominator:.4f}")
print(f"\nSo, P = 1 - {numerator:.4f} / {denominator:.4f}")
print(f"P = 1 - {prob_hit:.4f}")
print(f"P = {prob_never_hit:.4f}")

# Final result with two significant digits
print(f"\nThe approximate answer with two significant digits is: {prob_never_hit:.2f}")

# The required final output format
print(f"\n<<<{prob_never_hit:.2f}>>>")