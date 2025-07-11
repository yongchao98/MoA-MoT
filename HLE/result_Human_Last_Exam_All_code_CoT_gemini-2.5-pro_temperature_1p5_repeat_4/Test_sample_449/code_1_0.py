import numpy as np

# Step 1: Define constants and the starting point coordinates.
x, y = 3000, 4000
gamma = np.euler_gamma  # Euler-Mascheroni constant

# Step 2: Calculate the distance from the origin.
r = np.sqrt(x**2 + y**2)

# Step 3: Calculate the potential kernel a(x_0) for the starting point.
# The asymptotic formula for a(r) is (2/pi) * (ln(r) + gamma + (5/2)*ln(2)).
ln_r = np.log(r)
ln_2 = np.log(2)
# The constant part of the formula
K_term = gamma + (5/2) * ln_2
a_r = (2 / np.pi) * (ln_r + K_term)

# Step 4: The probability of the conditioned walk hitting the target set A is approximately 1/a(r).
# The probability of never hitting A is 1 - 1/a(r).
prob = 1 - 1 / a_r

# Step 5: Print the results following the required format.
# "Remember in the final code you still need to output each number in the final equation!"
# The equation is: 1 - 1 / a(r)
# where a(r) = (2/pi) * (ln(r) + gamma + (5/2)*ln(2))
print(f"The starting point is ({x},{y}), with distance r = {r:.0f} from the origin.")
print("The probability is calculated as: 1 - 1 / a(r)")
print(f"The value of the potential kernel a(r) is approximately {a_r:.4f}")
print("The calculation is:")
print(f"P(never hit A) = 1 - 1 / ({a_r:.4f})")
print(f"Result = {prob:.4f}")
print(f"The approximate answer with two significant digits is {prob:.2f}")

# Final answer block
final_answer = round(prob, 2)
# print(f"<<<{final_answer}>>>")