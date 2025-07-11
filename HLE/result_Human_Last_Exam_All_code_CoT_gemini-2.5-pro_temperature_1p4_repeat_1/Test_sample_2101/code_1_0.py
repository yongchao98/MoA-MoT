import math

# The problem is to find the difference between the probability of an electron
# escaping through the hypotenuse and the probability of it escaping through
# either of the two legs of an isosceles right triangle.

# The angles of an isosceles right triangle are 90, 45, and 45 degrees.
# In radians, these are pi/2, pi/4, and pi/4.
angle_opposite_hypotenuse = math.pi / 2.0
angle_opposite_leg = math.pi / 4.0

# The probability (P) of escape through a side is given by the formula:
# P = (pi - V) / (2 * pi)
# where V is the angle of the vertex opposite to that side.

# Calculate the probability of escaping through the hypotenuse.
# The angle opposite the hypotenuse is pi/2.
p_hyp = (math.pi - angle_opposite_hypotenuse) / (2 * math.pi)

# Calculate the probability of escaping through one of the legs.
# The angle opposite each leg is pi/4.
p_leg = (math.pi - angle_opposite_leg) / (2 * math.pi)

# The probability of escaping through either of the two legs is the sum,
# as they are mutually exclusive events.
p_legs_total = p_leg + p_leg

# Calculate the required difference.
difference = p_hyp - p_legs_total

# Print the results of the calculation step-by-step.
print("1. Probability of escaping through the hypotenuse (P_hyp):")
print(f"   P_hyp = (pi - pi/2) / (2*pi) = (pi/2) / (2*pi) = 1/4 = {p_hyp}")

print("\n2. Probability of escaping through one leg (P_leg):")
print(f"   P_leg = (pi - pi/4) / (2*pi) = (3*pi/4) / (2*pi) = 3/8 = {p_leg}")

print("\n3. Total probability of escaping through both legs (P_legs):")
print(f"   P_legs = P_leg + P_leg = 3/8 + 3/8 = 6/8 = 3/4 = {p_legs_total}")

print("\n4. The difference is P_hyp - P_legs:")
# The final equation with each number printed out
print(f"   {p_hyp} - {p_legs_total} = {difference}")

print("\nIn fraction form, the final equation is:")
print(f"   {1}/{4} - {3}/{4} = {-2}/{4} = {-1}/{2}")
print(f"\nThus, the final difference is {difference}.")

<<< -0.5 >>>