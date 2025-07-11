import numpy as np

def infinity_norm(v):
  """Computes the L-infinity norm of a vector."""
  return np.max(np.abs(v))

# 2. Define the endpoints of the line segment
a = np.array([1, 0])
b = np.array([1, 1])

# 3. Verify that a and b are on the unit sphere
norm_a = infinity_norm(a)
norm_b = infinity_norm(b)

print(f"Point a = {a}")
print(f"Point b = {b}")
print(f"Norm of a: ||a||_inf = {norm_a}")
print(f"Norm of b: ||b||_inf = {norm_b}")
print("Both points are on the unit sphere.")
print("-" * 20)

# 4. Verify that the line segment L(a,b) is on the unit sphere.
# A point on the segment is z(t) = (1-t)a + t*b
# For a=(1,0) and b=(1,1), z(t) = (1, t) for t in [0,1].
# We check a few sample points.
print("Verifying segment is on the unit sphere:")
for t in [0, 0.25, 0.5, 0.75, 1.0]:
  z = (1 - t) * a + t * b
  norm_z = infinity_norm(z)
  print(f"  For t={t:.2f}, point z = {z}, norm ||z||_inf = {norm_z}")
# For z=(1,t) with t in [0,1], max(|1|,|t|) is always 1.
print("The entire segment lies on the unit sphere.")
print("-" * 20)

# 5. Calculate the length of the segment, L = ||a-b||
segment_length = infinity_norm(a - b)

# Output the final equation with all numbers
a_minus_b = a-b
print("The length of the segment is ||a-b||_inf.")
print(f"a - b = {a} - {b} = {a_minus_b}")
print(f"||{a_minus_b}||_inf = {segment_length}")
print("\nThe final equation for the length is:")
print(f"|| {list(a)} - {list(b)} || = {segment_length}")
