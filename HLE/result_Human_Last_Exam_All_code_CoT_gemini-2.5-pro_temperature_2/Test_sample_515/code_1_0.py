import math

# Step 1: Define the constants and parameters
# The starting point is x0 = (0,1)
# The target point is z = (3600,0)
z_norm = 3600.0

# Mathematical constants
gamma = 0.57721566490153286060
# ln(2*sqrt(2)) = ln(sqrt(8)) = 0.5 * ln(8)
ln_2_sqrt_2 = 0.5 * math.log(8)

# Step 2: Calculate the potential kernel values
# a(x0) for x0 = (0,1) is known to be 4/pi
a_x0 = 4.0 / math.pi

# a(z) is calculated using the asymptotic formula for large ||z||
# a(z) = (2/pi) * (ln(||z||) + gamma + ln(2*sqrt(2)))
a_z = (2.0 / math.pi) * (math.log(z_norm) + gamma + ln_2_sqrt_2)

# Step 3: Calculate the probability
# The probability is the ratio a(x0) / a(z)
prob = a_x0 / a_z

# The formula simplifies to 2 / (ln(||z||) + gamma + ln(2*sqrt(2)))
denominator = math.log(z_norm) + gamma + ln_2_sqrt_2
prob_simplified = 2.0 / denominator

# Step 4: Print the results and the final equation
print("The problem is to find the probability that a random walk, conditioned to avoid the origin,")
print("starting at (0,1), hits the set of neighbors of z = (3600,0).")
print("\nThis probability is well-approximated by P(hit z before origin), which is given by the ratio of potential kernels a(x0)/a(z).")
print("\nCalculation steps:")
print(f"1. Potential kernel at starting point x0=(0,1): a(x0) = 4/\u03C0 \u2248 {a_x0:.4f}")
print(f"2. Potential kernel at target point z=(3600,0) using asymptotic formula:")
print(f"   a(z) \u2248 (2/\u03C0) * (ln({z_norm}) + \u03B3 + ln(2\u221A2))")
print(f"   ln({z_norm}) \u2248 {math.log(z_norm):.4f}")
print(f"   \u03B3 \u2248 {gamma:.4f}")
print(f"   ln(2\u221A2) \u2248 {ln_2_sqrt_2:.4f}")
print(f"   a(z) \u2248 {a_z:.4f}")
print("\n3. The final probability is the ratio a(x0)/a(z):")
print(f"   P \u2248 a(x0)/a(z) = (4/\u03C0) / [(2/\u03C0) * (ln({z_norm}) + \u03B3 + ln(2\u221A2))]")
print(f"   P = 2 / (ln({z_norm}) + {gamma:.4f} + {ln_2_sqrt_2:.4f})")
print(f"   P = 2 / ({math.log(z_norm):.4f} + {gamma:.4f} + {ln_2_sqrt_2:.4f})")
print(f"   P = 2 / {denominator:.4f}")
print(f"   P \u2248 {prob:.4f}")
print(f"\nApproximate answer with two significant digits: {prob:.2g}")
