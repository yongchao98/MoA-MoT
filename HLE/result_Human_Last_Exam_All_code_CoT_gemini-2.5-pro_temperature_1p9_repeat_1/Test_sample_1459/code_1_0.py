import math

# Step 1: Define the parameters of the spaces
# H is the length of the interval [0, 1]
H = 1.0
# L is the circumference of the unit circle (radius=1)
L = 2 * math.pi

# Step 2: Verify the condition for the formula
# The formula d_GH = (L - 2H) / (2 * sqrt(2)) is valid for H/L <= 1/4.
ratio = H / L
condition_met = ratio <= 0.25

# Step 3: Define the numbers in the simplified formula (pi - 1) / sqrt(2)
pi_val = math.pi
one_val = 1.0
sqrt2_val = math.sqrt(2)

# Step 4: Calculate the numerator and the denominator
numerator = pi_val - one_val
denominator = sqrt2_val

# Step 5: Calculate the final Gromov-Hausdorff distance
gh_distance = numerator / denominator

# Step 6: Print the result showing each number in the equation
print(f"The Gromov-Hausdorff distance between the interval [0,1] and the unit circle is given by the formula (pi - 1) / sqrt(2).")
print("This is calculated as:")
print(f"({pi_val} - {one_val}) / {sqrt2_val} = {gh_distance}")
