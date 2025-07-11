import math

# Step 1: Define constants and the starting point
x0 = (3000, 4000)
a1 = 1.0  # a(1,0) with normalization a(0,0)=0

# Step 2: Calculate the distance of the starting point from the origin
r = math.sqrt(x0[0]**2 + x0[1]**2)

# Step 3: Calculate the constant C in the asymptotic formula for a(r)
# C = (2/pi) * (gamma_E + ln(8))
gamma_E = 0.57721566490153287
C = (2 / math.pi) * (gamma_E + math.log(8))

# Step 4: Calculate a(x0) using the asymptotic formula
# a(r) ~ (2/pi) * ln(r) + C
a_x0 = (2 / math.pi) * math.log(r) + C

# Step 5: Calculate the final probability
# P = 1 - a(1,0) / a(x0)
prob = 1 - a1 / a_x0

# Print the final result and the intermediate steps for clarity
print(f"The calculation is for the probability P = 1 - a(1,0) / a({x0[0]},{x0[1]}).")
print(f"The value of a(1,0) is {a1:.1f}.")
print(f"The distance of the starting point from the origin is |x0| = {r:.0f}.")
print(f"The value of a({x0[0]},{x0[1]}) is approximated by a({r:.0f}) ≈ {a_x0:.4f}.")
print(f"So the equation is: P = 1 - {a1:.1f} / {a_x0:.4f}")
print(f"The final probability is P ≈ {prob:.4f}.")

# Provide the approximate answer with two significant digits.
# The probability 0.8594 becomes 0.86.
final_answer = round(prob, 2)
print(f"The approximate answer with two significant digits is: {final_answer}")
<<<0.86>>>