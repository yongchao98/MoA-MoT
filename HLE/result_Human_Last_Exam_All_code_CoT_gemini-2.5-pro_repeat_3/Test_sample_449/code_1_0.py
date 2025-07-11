import math

# Step 1: Define constants and the starting point's norm
euler_mascheroni = 0.5772156649
pi = math.pi
x0_norm = 5000.0

# Step 2: Calculate the value of the potential kernel at the origin, a_plus(0)
# a_plus(0) = (2/pi) * (gamma + 2*ln(2))
a_plus_0 = (2 / pi) * (euler_mascheroni + 2 * math.log(2))

# Step 3: Calculate the value of the potential kernel at a neighbor of the origin, e.g., (1,0)
# a_plus(1,0) = a_plus(0) - 1
a_plus_1 = a_plus_0 - 1

# Step 4: Calculate the value of the potential kernel at the starting point (3000,4000)
# a_plus(x0) is approximated by (2/pi)*ln(||x0||) + a_plus(0)
a_plus_x0 = (2 / pi) * math.log(x0_norm) + a_plus_0

# Step 5: Calculate the final probability
# P = 1 - a_plus(1,0) / a_plus(x0)
probability = 1 - (a_plus_1 / a_plus_x0)

# Output the results in a clear format
print(f"The potential kernel at the origin, a^+(0), is approximately: {a_plus_0:.4f}")
print(f"The potential kernel at a neighbor of the origin, a^+(1,0), is approximately: {a_plus_1:.4f}")
print(f"The potential kernel at the starting point (3000,4000) is approximately: {a_plus_x0:.4f}")
print(f"The probability is calculated as: 1 - ({a_plus_1:.4f} / {a_plus_x0:.4f})")
print(f"The final probability is approximately: {probability:.4f}")

# Approximate answer with two significant digits
# The result 0.9625... rounds to 0.96
final_answer = round(probability, 2)
print(f"\nThe approximate answer with two significant digits is: {final_answer}")
<<<0.96>>>