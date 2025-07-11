import math

# Step 1: Define the parameters of the two metric spaces.
# Space X: The interval [0, 1]. Its length is L.
L = 1.0

# Space Y: The unit circle, interpreted as a circle of circumference 1.
C = 1.0

print(f"Calculating the Gromov-Hausdorff distance between:")
print(f"- An interval of length L = {L}")
print(f"- A circle of circumference C = {C}\n")

# Step 2: Check the condition for the formula to be used.
# The simpler formula d_GH = L/2 applies if L >= C/pi.
C_div_pi = C / math.pi
condition_is_met = L >= C_div_pi

print(f"The condition for the simpler formula is L >= C/π.")
print(f"Checking the condition: {L} >= {C}/π ≈ {C_div_pi:.4f}")
print(f"Is the condition met? {condition_is_met}\n")


# Step 3: Apply the correct formula and calculate the distance.
if condition_is_met:
    # Use the formula d_GH = L/2
    distance = L / 2
    print(f"The condition is met. The applicable formula is d_GH = L / 2.")
    print(f"Final Equation: {L} / 2 = {distance}")
else:
    # This case is not reached for these specific spaces, but included for completeness.
    # The formula for L < C/pi is more complex and not needed here.
    print("The condition for the simpler formula was not met.")
    distance = "Calculation requires a more complex formula."
    print(f"Result: {distance}")

print(f"\nThe Gromov-Hausdorff distance is: {distance}")
