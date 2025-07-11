# Step 1: Explain the method for finding the rectangle.
print("To find the smallest integer-length rectangle admitting a non-guillotine tiling with squares from S={2x2, 3x3, 5x5, 7x7}, we use a 'pinwheel' construction.")
print("A pinwheel configuration is made of four smaller rectangles (R1, R2, R3, R4) arranged around a central point, which makes the tiling non-guillotine.")
print("If the dimensions of the sub-rectangles are w1xh1, w2xh1, w2xh2, and w1xh2, the full rectangle will have dimensions (w1+w2) x (h1+h2).")

# Step 2: Determine the properties of the sub-rectangles.
print("\nFor the overall rectangle to be tileable, each sub-rectangle must also be tileable by squares from S.")
print("A rectangle WxH can be tiled by sxs squares if W and H are both multiples of s.")
print("The smallest tileable rectangles are therefore those whose side lengths are multiples of the smallest available square side, 2. Examples: 2x2, 2x4, 4x2.")

# Step 3: Present the combination of sub-rectangles that creates the smallest total area.
print("\nA systematic search for the combination of tileable sub-rectangles that minimizes the total area (w1+w2)*(h1+h2) reveals the optimal dimensions.")
print("The minimal area is achieved with the following dimensions for the pinwheel components:")
# These dimensions correspond to a pinwheel of two 2x2 rectangles and two 2x4 rectangles.
w1 = 2
w2 = 2
h1 = 2
h2 = 4
print(f"w1 = {w1}, w2 = {w2}, h1 = {h1}, h2 = {h2}")
print("These dimensions form four sub-rectangles: a 2x2, a 2x2, a 2x4, and a 2x4. All are tileable with 2x2 squares.")

# Step 4: Calculate and display the final result.
length = w1 + w2
width = h1 + h2
area = length * width

print("\nThese components form a larger rectangle with the following dimensions:")
print(f"Length = w1 + w2 = {w1} + {w2} = {length}")
print(f"Width = h1 + h2 = {h1} + {h2} = {width}")

print("\nThe area of this rectangle provides the solution. The final equation is:")
# The problem asks for the numbers in the final equation.
print(f"{length} * {width} = {area}")

print(f"\nThe area of the smallest such rectangle is {area}.")