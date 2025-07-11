import math

# This script calculates the difference between the probability that an electron
# escapes through the hypotenuse and the probability that it escapes through
# either of the two legs of an isosceles right triangle.

# As outlined in the plan, the problem boils down to calculating the average
# angle subtended by the hypotenuse at a random point inside the triangle.

# Step 1: Use the known analytical result for the average angle subtended by the
# hypotenuse in an isosceles right triangle.
# The result is: <angle_h> = pi/2 + ln(2) - 1
avg_angle_h = math.pi / 2 + math.log(2) - 1

# Step 2: Calculate the probability of escaping through the hypotenuse (P_h).
# This is the average angle divided by the total angle of 2*pi.
P_h = avg_angle_h / (2 * math.pi)

# Step 3: The probability of escaping through the legs (P_legs) is the rest of the
# total probability.
P_legs = 1 - P_h

# Step 4: The final difference is D = P_h - P_legs.
difference = P_h - P_legs

# Output the results of the calculation step-by-step.
# The final equation is Difference = P_h - P_legs
print("The probability of escaping through the hypotenuse is:")
print(f"P_h = {P_h:.4f}")
print("\nThe probability of escaping through the legs is:")
print(f"P_legs = {P_legs:.4f}")
print("\nThe difference is calculated as:")
print(f"Difference = {P_h:.4f} - {P_legs:.4f} = {difference:.4f}")