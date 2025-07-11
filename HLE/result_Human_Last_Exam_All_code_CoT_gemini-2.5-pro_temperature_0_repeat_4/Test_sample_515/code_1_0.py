import math

# Step 1: Define the constants
# Euler-Mascheroni constant
gamma = 0.5772156649
# Pi
pi = math.pi
# Natural logarithm of 2
ln2 = math.log(2)

# Step 2: Define the coordinates
# Starting point A = (0,1)
# Target point B = (3600,0)
x_B = 3600
y_B = 0

# Step 3: Calculate the potential kernel a(A) and a(B)
# a(A) = a(0,1) is known to be exactly 1.
a_A = 1.0
print(f"The potential kernel at the starting point A=(0,1) is a(A) = {a_A}")

# Calculate the norm of B
norm_B = math.sqrt(x_B**2 + y_B**2)

# Calculate a(B) using the asymptotic formula
# a(B) ~ (2/pi) * (ln|B| + gamma + 3*ln2)
ln_norm_B = math.log(norm_B)
a_B = (2 / pi) * (ln_norm_B + gamma + 3 * ln2)

print(f"The norm of the target point B=(3600,0) is |B| = {norm_B}")
print(f"The natural logarithm of the norm is ln|B| = {ln_norm_B:.4f}")
print(f"The constant part of the formula is gamma + 3*ln(2) = {gamma + 3*ln2:.4f}")
print(f"The potential kernel at the target point B is approximately a(B) = {a_B:.4f}")

# Step 4: Calculate the probability
# The probability is approximated by the ratio a(A) / a(B)
prob = a_A / a_B

print(f"\nThe probability is p = a(A) / a(B) = {a_A} / {a_B:.4f} = {prob:.4f}")

# Step 5: Round to two significant digits
# To round to two significant digits, we can use a format string.
prob_rounded = float(f"{prob:.2g}")
print(f"The approximate answer with two significant digits is {prob_rounded}")
