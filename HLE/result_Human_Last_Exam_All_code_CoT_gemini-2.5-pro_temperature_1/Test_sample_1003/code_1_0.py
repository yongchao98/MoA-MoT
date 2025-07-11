import math

# Step 1: Define the known values from the second frame of reference.
# Angle between S1 and S2 is pi/2.
cos_theta_12 = 0
# Angle between S2 and S3 is 3*pi/4.
cos_theta_23 = -1 / math.sqrt(2)

# Step 2: The desired ratio is given by the formula derived from the Lorentz invariant cross-ratio:
# (1 - cos(theta_14)) / (1 - cos(theta_34)) = (1 - cos(theta_12)) / (1 - cos(theta_23))
# We will now substitute the known values into the right side of the equation.

numerator = 1 - cos_theta_12
denominator = 1 - cos_theta_23

# Step 3: Print the equation with the numbers plugged in.
print("The final equation is:")
print(f"(1 - cos(theta_14)) / (1 - cos(theta_34)) = (1 - {cos_theta_12}) / (1 - ({cos_theta_23:.4f}))")
print(f"                             = {numerator} / {denominator:.4f}")

# Step 4: Calculate the final value.
result = numerator / denominator

print(f"\nThe numerical value of the ratio is: {result:.4f}")
print(f"The exact value is 2 - sqrt(2).")

# Final answer in the specified format
final_answer = 2 - math.sqrt(2)
# The problem statement does not specify a format for the final answer, so we will provide the numerical value.
# <<<2 - math.sqrt(2)>>> is not a valid format.
# We will output the numerical value.
# Let's re-read the instructions. "return the answer with the format <<<answer content>>> at the end of your response".
# The content can be a number or a letter. So a float is fine.