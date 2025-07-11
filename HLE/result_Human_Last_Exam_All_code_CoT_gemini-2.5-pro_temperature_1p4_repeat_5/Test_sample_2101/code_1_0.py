import math

# Step 1: Define the relative lengths of the triangle's sides.
# We can set the leg length L=1.0 without loss of generality, as it will
# cancel out when calculating probabilities.
leg_length = 1.0
hypotenuse_length = 1.0 * math.sqrt(2)

# Step 2: Calculate the total perimeter.
perimeter = leg_length + leg_length + hypotenuse_length

# Step 3: Use the principle from integral geometry to find the probabilities.
# P(escape through a side) = (length of the side) / (total perimeter)
prob_hypotenuse = hypotenuse_length / perimeter
prob_leg = leg_length / perimeter

# Step 4: Calculate the required difference.
# D = P(hypotenuse) - (P(leg1) + P(leg2))
difference = prob_hypotenuse - (prob_leg + prob_leg)

# Step 5: Print the final equation with the calculated numbers.
# The problem asks for the difference between the probability of escaping through the
# hypotenuse and the probability of escaping through either of the two legs.
print("Let P_hyp be the probability of escaping through the hypotenuse.")
print("Let P_leg be the probability of escaping through one leg.")
print("We need to find D = P_hyp - (P_leg + P_leg).")
print("\nThe values are:")
print(f"P_hyp = {prob_hypotenuse}")
print(f"P_leg = {prob_leg}")
print("\nThe final equation with these numbers is:")
print(f"{prob_hypotenuse} - ({prob_leg} + {prob_leg}) = {difference}")
print("\nThe exact analytical result is 2*sqrt(2) - 3.")
