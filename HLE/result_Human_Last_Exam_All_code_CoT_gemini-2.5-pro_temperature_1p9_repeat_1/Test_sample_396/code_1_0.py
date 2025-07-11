import numpy as np

# Step 1 & 2: Decrypt the message and define the initial vectors.
# "nggyu..." -> UP -> (0, 0, 1)
# "...nglyd..." -> DOWN -> (0, 0, -1)
# "...ngraady..." -> RIGHT (from 'r' in 'run') -> (1, 0, 0)
# Z is the vertical axis.
initial_vectors = [
    (0, 0, 1),  # UP
    (0, 0, -1), # DOWN
    (1, 0, 0)   # RIGHT
]

final_vectors = []

# Step 3: Iterate through vectors and apply rotation where needed.
# We process a 1-based list, so the 2nd vector has an index of 1.
for i, vec in enumerate(initial_vectors):
    # Check if it's the 2nd, 4th, etc., vector.
    if (i + 1) % 2 == 0:
        x, y, z = vec
        # Apply a 90-degree clockwise rotation around the x-axis.
        # The transformation is (x, y, z) -> (x, -z, y).
        rotated_vec = (x, -z, y)
        final_vectors.append(rotated_vec)
    else:
        final_vectors.append(vec)

# Step 4: Sum the final vectors.
# We use numpy for a concise sum, but a manual sum would also work.
sum_vector = np.sum(final_vectors, axis=0)
sum_x, sum_y, sum_z = sum_vector

# Step 5: Format the output to show the summation equation.
# The `str(v[c]).replace("-", "- ")` part ensures spacing for negative numbers.
x_components = " + ".join(str(v[0]) for v in final_vectors)
y_components = " + ".join(str(v[1]) for v in final_vectors)
z_components = " + ".join(str(v[2]) for v in final_vectors)

# Print the final result in the requested format.
print(f"({x_components}, {y_components}, {z_components}) = ({sum_x},{sum_y},{sum_z})")