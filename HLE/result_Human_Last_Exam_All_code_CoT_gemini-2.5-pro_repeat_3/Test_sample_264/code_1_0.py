# Step 1: Define the properties of the cubes and the shift.
# We place the first cube (C1) with one corner at the origin.
c1_min = 0.0
c1_max = 1.0

# The shift is 1/2 along the main diagonal, which we interpret as a vector (0.5, 0.5, 0.5).
shift = 0.5

# Step 2: Define the boundaries of the second cube (C2) after the shift.
c2_min = c1_min + shift
c2_max = c1_max + shift

# Step 3: Find the boundaries of the intersection volume.
# The intersection is a cube itself. We calculate the start and end of its sides.
inter_min = max(c1_min, c2_min)
inter_max = min(c1_max, c2_max)

# Step 4: Calculate the side length of the intersection cube.
side_length = inter_max - inter_min

# Step 5: Calculate the volume of the intersection.
volume = side_length ** 3

print("This script calculates the volume of the intersection of two cubes.")
print("Each cube has a side length of 1.")
print("One cube is shifted by 1/2 along the main diagonal relative to the other.")
print("\nThe intersection forms a smaller cube.")
print(f"The side length of this intersection cube is: {inter_max} - {inter_min} = {side_length}")
print("\nThe volume is the cube of its side length.")
# The final equation is printed with each number, as requested.
print(f"Volume = {side_length} * {side_length} * {side_length} = {volume}")
print("\nIn fraction form:")
print("Volume = (1/2) * (1/2) * (1/2) = 1/8")
