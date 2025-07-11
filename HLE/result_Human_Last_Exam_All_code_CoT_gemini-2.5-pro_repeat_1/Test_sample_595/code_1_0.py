import math

# Define the lengths of the triangle sides
leg_length = 18

# Choose the x-component for the first leg vector
x1 = 15
# Calculate the y-component
y1 = math.sqrt(leg_length**2 - x1**2)

# Place the right-angle vertex A at a non-lattice point
A = (0.1, 0.1)

# Calculate the coordinates of the other vertices B and C
# B = A + (x1, y1)
# C = A + (-y1, x1)
B = (A[0] + x1, A[1] + y1)
C = (A[0] - y1, A[1] + x1)

# Function to calculate the number of squares crossed by a segment
def squares_crossed(p1, p2):
    x_min, x_max = min(p1[0], p2[0]), max(p1[0], p2[0])
    y_min, y_max = min(p1[1], p2[1]), max(p1[1], p2[1])
    
    # Number of vertical grid lines crossed
    # This is the number of integers in the open interval (x_min, x_max)
    vert_crossings = math.floor(x_max) - math.ceil(x_min) + 1
    # Adjust for endpoints being very close to integers
    if x_max == math.floor(x_max): vert_crossings -=1
    if x_min == math.ceil(x_min): vert_crossings -=1

    # Number of horizontal grid lines crossed
    horiz_crossings = math.floor(y_max) - math.ceil(y_min) + 1
    if y_max == math.floor(y_max): horiz_crossings -=1
    if y_min == math.ceil(y_min): horiz_crossings -=1
    
    return vert_crossings + horiz_crossings

# Calculate crossings for each side
k_AB = squares_crossed(A, B)
k_AC = squares_crossed(A, C)
k_BC = squares_crossed(B, C)

# Calculate the total
k_total = k_AB + k_AC + k_BC

print("Calculation Steps:")
print(f"Chosen vector components for one leg: ({x1}, sqrt({leg_length**2 - x1**2})) = ({x1}, {y1:.2f})")
print(f"Vertex A: {A}")
print(f"Vertex B: ({B[0]:.2f}, {B[1]:.2f})")
print(f"Vertex C: ({C[0]:.2f}, {C[1]:.2f})")
print("-" * 20)
print(f"Squares crossed by leg AB: {k_AB}")
print(f"Squares crossed by leg AC: {k_AC}")
print(f"Squares crossed by hypotenuse BC: {k_BC}")
print("-" * 20)
print("The final equation is the sum of the squares crossed by each side:")
print(f"k = {k_AB} + {k_AC} + {k_BC}")
print(f"k = {k_total}")
