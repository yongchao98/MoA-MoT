import math

# An isosceles right triangle has angles of 90, 45, and 45 degrees.
# In radians, these are pi/2, pi/4, and pi/4.
angle_opposite_hypotenuse = math.pi / 2
angle_opposite_leg = math.pi / 4

# A theorem in geometric probability states that the probability of a random ray
# from a random point inside a triangle escaping through a given side is
# the angle of the opposite vertex divided by pi.

# Calculate the probability of the electron escaping through the hypotenuse.
# The hypotenuse is opposite the 90-degree (pi/2 radians) angle.
p_hyp = angle_opposite_hypotenuse / math.pi

# Calculate the probability of the electron escaping through the two legs.
# Each leg is opposite a 45-degree (pi/4 radians) angle.
p_leg_1 = angle_opposite_leg / math.pi
p_leg_2 = angle_opposite_leg / math.pi
# The total probability is the sum for both legs.
p_legs = p_leg_1 + p_leg_2

# Calculate the difference between the two probabilities.
difference = p_hyp - p_legs

# Print the results step-by-step
print(f"The probability of escaping through the hypotenuse is calculated as the opposite angle (pi/2) divided by pi.")
print(f"P(hypotenuse) = (pi/2) / pi = {p_hyp}")
print("")
print(f"The probability of escaping through one leg is the opposite angle (pi/4) divided by pi.")
print(f"P(one leg) = (pi/4) / pi = {p_leg_1}")
print(f"The total probability of escaping through either leg is the sum of the probabilities for each leg.")
print(f"P(legs) = {p_leg_1} + {p_leg_2} = {p_legs}")
print("")
print(f"The difference between the probability of escaping through the hypotenuse and through the legs is:")
print(f"Difference = P(hypotenuse) - P(legs)")
print(f"Difference = {p_hyp} - {p_legs} = {difference}")