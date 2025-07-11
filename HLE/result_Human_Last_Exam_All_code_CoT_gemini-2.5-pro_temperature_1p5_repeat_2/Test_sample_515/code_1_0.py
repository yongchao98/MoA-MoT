import math

# Parameters for the calculation
# Starting point x0 = (0,1), so |x0| = 1
# Target point z = (3600,0), so |z| = 3600

# Value of the potential kernel at x0=(0,1)
# a_x0 = 2 / math.pi
# We can see from the formula that a_x0 will cancel out with the 2/pi factor in a_z.
# The probability is 1 / (log|z| + C') where C' is the constant part.

z_norm = 3600.0
euler_gamma = 0.5772156649

# The probability is given by the formula 1 / (log(|z|) + gamma + 2*log(2))
log_z = math.log(z_norm)
log_2 = math.log(2)

# The constant part of the denominator
constant_part = euler_gamma + 2 * log_2

# Denominator of the probability expression
denominator = log_z + constant_part

# The final probability
probability = 1 / denominator

print(f"The position z is ({int(z_norm)}, 0)")
print(f"The potential kernel at x0=(0,1) is a(x0) = 2/pi")
print(f"The potential kernel at z is approximately a(z) = (2/pi) * (log({int(z_norm)}) + gamma + 2*log(2))")
print(f"The probability is P = a(x0) / a(z) = 1 / (log({int(z_norm)}) + gamma + 2*log(2))")
print(f"log({int(z_norm)}) = {log_z:.4f}")
print(f"gamma = {euler_gamma:.4f}")
print(f"2*log(2) = {2*log_2:.4f}")
print(f"P = 1 / ({log_z:.4f} + {euler_gamma:.4f} + {2*log_2:.4f}) = 1 / {denominator:.4f}")
print(f"P = {probability:.4f}")

# Approximate answer with two significant digits
approx_prob = round(probability, 2)
# Correction for a number like 0.098 which has two sig figs.
if probability < 0.1 and probability > 0.01:
    approx_prob_str = f"{probability:.2g}"
else:
    approx_prob_str = f"{probability:.2f}"

print(f"\nThe approximate answer with two significant digits is: {approx_prob_str}")