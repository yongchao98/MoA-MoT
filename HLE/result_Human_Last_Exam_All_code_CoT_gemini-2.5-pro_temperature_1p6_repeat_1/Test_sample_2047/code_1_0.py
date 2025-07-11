import math

# I will use d=4 and lambda=1 as they are the smallest integer values
# for which the function is defined as per the problem description (d>=4, lambda>=1).
d = 4
lambda_val = 1

# The derived formula for l(d, lambda) is:
# l(d, lambda) = (1 / (2 * lambda)) * [arccos(sqrt(2/d))^2 - arccos(sqrt(3/d))^2]

# First, calculate the arguments for the arccos function
# For the second term in the expression: arccos(sqrt(3/d))
val_sqrt_3_d = math.sqrt(3 / d)
# For the first term in the expression: arccos(sqrt(2/d))
val_sqrt_2_d = math.sqrt(2 / d)

# Second, calculate the angles (theta1 and theta2)
# theta1 corresponds to x1
theta1 = math.acos(val_sqrt_3_d)
# theta2 corresponds to x2
theta2 = math.acos(val_sqrt_2_d)

# Third, square these angles
theta1_sq = theta1**2
theta2_sq = theta2**2

# Finally, apply the full formula
result = (1 / (2 * lambda_val)) * (theta2_sq - theta1_sq)

print(f"Based on the problem analysis, the expression for l(d, lambda) is:")
print(f"l(d, lambda) = (1/(2*lambda)) * [arccos(sqrt(2/d))^2 - arccos(sqrt(3/d))^2]")
print(f"We will compute this for d = {d} and lambda = {lambda_val}.\n")

print("The final equation is: l(4, 1) = (1/(2*1)) * [arccos(sqrt(2/4))^2 - arccos(sqrt(3/4))^2]")
print("Let's substitute the numbers step-by-step:\n")

print(f"For the arccos(sqrt(3/4)) term:")
print(f"The number sqrt(3/4) is: {val_sqrt_3_d:.8f}")
print(f"The number arccos({val_sqrt_3_d:.8f}) is: {theta1:.8f} radians (which is pi/6)")
print(f"The number (arccos(sqrt(3/4)))^2 is: {theta1_sq:.8f}\n")

print(f"For the arccos(sqrt(2/4)) term:")
print(f"The number sqrt(2/4) is: {val_sqrt_2_d:.8f}")
print(f"The number arccos({val_sqrt_2_d:.8f}) is: {theta2:.8f} radians (which is pi/4)")
print(f"The number (arccos(sqrt(2/4)))^2 is: {theta2_sq:.8f}\n")

print("Putting it all together:")
print(f"l(4, 1) = (1/2) * [{theta2_sq:.8f} - {theta1_sq:.8f}]")
print(f"l(4, 1) = 0.5 * [{theta2_sq - theta1_sq:.8f}]")
print(f"The final value is: {result:.8f}")
