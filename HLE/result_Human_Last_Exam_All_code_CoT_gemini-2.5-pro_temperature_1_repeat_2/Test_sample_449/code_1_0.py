import math

# Step 1: Define the constants and the starting point
x0 = (3000, 4000)
gamma = 0.5772156649  # Euler-Mascheroni constant

# Step 2: Calculate the distance of the starting point from the origin
r0 = math.sqrt(x0[0]**2 + x0[1]**2)

# Step 3: Calculate the values of the potential kernel a(x)
# For n = (1,0), the exact value is a(1,0) = 2/pi
a_n = 2 / math.pi

# For the starting point x0, use the asymptotic formula for large |x|:
# a(x) ~ (1/pi) * (ln|x| + gamma + ln(8))
a_x0 = (1 / math.pi) * (math.log(r0) + gamma + math.log(8))

# Step 4: Calculate the probability q(x0) of hitting the set N
# q(x0) = a(n) / a(x0)
q_x0 = a_n / a_x0

# Step 5: The desired probability is the probability of never hitting N, which is 1 - q(x0)
prob_never_hit_N = 1 - q_x0

# Print the final equation with all the numbers
print(f"The probability is given by the expression: 1 - a(1,0) / a(3000,4000)")
print(f"a(1,0) = 2 / pi = {a_n:.4f}")
print(f"a(3000,4000) ~= (1/pi) * (ln({r0}) + gamma + ln(8)) = {a_x0:.4f}")
print(f"Probability = 1 - {a_n:.4f} / {a_x0:.4f}")
print(f"Probability = 1 - {q_x0:.4f}")
print(f"Probability = {prob_never_hit_N:.4f}")

# Print the final answer with two significant digits
print(f"\nThe approximate answer with two significant digits is: {prob_never_hit_N:.2f}")
