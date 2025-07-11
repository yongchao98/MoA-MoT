import math

# Step 1: Define the parameters based on the problem description.
N = 5  # Number of black keys (self-similar copies)

# Dimensions of the whole keyboard
W_total = 3.0
H_total = 1.0

# Dimensions of the keys
num_white_keys = 7
W_white = W_total / num_white_keys
W_black = 0.5 * W_white
H_black = (9.0 / 14.0) * H_total

# Step 2: Calculate the scaling factors.
# s_x is the ratio of the new width (a black key) to the original width.
s_x = W_black / W_total
# s_y is the ratio of the new height (a black key) to the original height.
s_y = H_black / H_total

# Step 3: Check the condition for the dimension formula D = D_x + 1.
# The dimension of the projection on the x-axis is D_x = log(N) / log(1/s_x).
# The condition simplifies to s_y >= s_x.
condition_met = s_y >= s_x

# Step 4: Calculate the Minkowski-Bouligand dimension.
# D = D_x + 1 = (log(N) / log(1/s_x)) + 1
inv_sx = 1.0 / s_x
Dx = math.log(N) / math.log(inv_sx)
D = Dx + 1

# Step 5: Print the explanation and the final result.
print("Calculating the Minkowski-Bouligand Dimension (D) for the fractal piano keys.")
print("The dimension D for this self-affine fractal is given by the formula D = Dx + 1, where Dx is the dimension of the projection onto the x-axis.")
print(f"The number of recursive copies, N, is the number of black keys: {N}")
print(f"The scaling factor in the x-direction, s_x, is {s_x:.4f}")
print(f"The scaling factor in the y-direction, s_y, is {s_y:.4f}")
print(f"The condition for the formula (s_y >= s_x) is {'met' if condition_met else 'not met'}.")
print("\nThe dimension Dx is calculated as log(N) / log(1/s_x):")
print(f"Dx = log({N}) / log({inv_sx:.2f})")
print(f"Dx = {Dx:.4f}")
print("\nThe final dimension D is Dx + 1:")
print(f"D = {Dx:.4f} + 1")
print(f"D = {D:.4f}")
