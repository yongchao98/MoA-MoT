import math

# Step 1: Define the optimal width (W) and height (H) of the triangle's bounding box.
# These values are derived from optimizing the orientation of the triangle.
# W = 18 * (cos(t) + sin(t)) and H = 18 * cos(t) where tan(t) = 1/2.
# This gives cos(t) = 2/sqrt(5) and sin(t) = 1/sqrt(5).
# W = 18 * (2/sqrt(5) + 1/sqrt(5)) = 54/sqrt(5)
# H = 18 * (2/sqrt(5)) = 36/sqrt(5)

W_val = 54 / math.sqrt(5)
H_val = 36 / math.sqrt(5)

# Step 2: Calculate the number of integer coordinates spanned by W and H.
# Since W and H are irrational, we can position the triangle to make the number
# of spanned integers equal to ceil(W) and ceil(H).
Nv = math.ceil(W_val)
Nh = math.ceil(H_val)

# Step 3: Calculate the total number of squares crossed, k.
# k = 2 * (number of vertical lines) + 2 * (number of horizontal lines)
k = 2 * Nv + 2 * Nh

# Step 4: Print the final equation with all the numbers.
print(f"The optimal bounding box has width W = 54/√5 ≈ {W_val:.4f} and height H = 36/√5 ≈ {H_val:.4f}.")
print("The largest number of squares k is calculated as:")
print(f"k = 2 * (ceil(W) + ceil(H))")
print(f"k = 2 * (ceil({W_val:.4f}) + ceil({H_val:.4f}))")
print(f"k = 2 * ({Nv} + {Nh})")
print(f"k = 2 * ({Nv + Nh})")
print(f"k = {k}")
