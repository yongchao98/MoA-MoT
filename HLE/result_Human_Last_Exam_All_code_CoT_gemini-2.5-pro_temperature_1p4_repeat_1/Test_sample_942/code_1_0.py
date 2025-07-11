import numpy as np

def l_inf_norm(v):
  """Calculates the L-infinity norm of a vector."""
  # The L-infinity norm is the maximum absolute value of the vector's elements.
  return np.max(np.abs(v))

print("Step 1: Establishing a theoretical maximum length.")
print("Let a and b be the endpoints of a line segment on the unit sphere.")
print("This means ||a|| = 1 and ||b|| = 1.")
print("The length of the segment is L = ||a - b||.")
print("By the triangle inequality, L <= ||a|| + ||b|| = 1 + 1 = 2.")
print("So, the maximum possible length is 2.")
print("\n------------------------------------------------\n")

print("Step 2: Finding a space where this maximum is achieved.")
print("We will use R^2 with the L-infinity norm: ||(x,y)|| = max(|x|, |y|).")
print("In this space, the unit ball is the square [-1,1] x [-1,1], which is metrically convex.")
print("The unit sphere is the boundary of this square.")
print("\n------------------------------------------------\n")

print("Step 3: Demonstrating a segment of length 2 on the unit sphere.")
# Define the endpoints of a line segment on the boundary of the square.
# This segment is one of the sides of the square.
a = np.array([1.0, -1.0])
b = np.array([1.0, 1.0])

# Verify endpoints are on the unit sphere
norm_a = l_inf_norm(a)
norm_b = l_inf_norm(b)

print(f"Let point a = {a}. Its norm is ||a||_inf = max(|{a[0]}|, |{a[1]}|) = {norm_a}.")
print(f"Let point b = {b}. Its norm is ||b||_inf = max(|{b[0]}|, |{b[1]}|) = {norm_b}.")
print("Both points are on the unit sphere.")
print()

# A point on the line segment is x(t) = (1-t)a + tb.
# For a=(1,-1) and b=(1,1), x(t) = (1, 2t-1) for t in [0,1].
# The norm is ||x(t)||_inf = max(|1|, |2t-1|). Since |2t-1|<=1, the norm is always 1.
# This confirms the segment from a to b is on the unit sphere.
print("The line segment connecting a and b also lies on the unit sphere.")
print("\n------------------------------------------------\n")

print("Step 4: Calculating the length of the segment.")
print("The length is ||a - b||_inf.")

# Perform the calculation
segment_vector = a - b

# Print the final equation with numbers as requested
print(f"Length = ||{a} - {b}||_inf")
print(f"       = ||({a[0] - b[0]}, {a[1] - b[1]})||_inf")
print(f"       = ||{segment_vector}||_inf")
print(f"       = max(|{segment_vector[0]}|, |{segment_vector[1]}|)")
print(f"       = {l_inf_norm(segment_vector)}")

print("\nSince the maximum possible length is 2, and we found an example of length 2, the answer is 2.")
