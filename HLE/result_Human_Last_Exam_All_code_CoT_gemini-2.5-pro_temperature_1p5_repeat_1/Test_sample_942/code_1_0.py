import numpy as np

def l_inf_norm(v):
    """Calculates the L-infinity norm of a vector."""
    if not isinstance(v, np.ndarray):
        v = np.array(v)
    return np.max(np.abs(v))

print("Demonstrating the case where the segment length is 2.")
print("We use the vector space R^2 with the L-infinity norm.")
print("The unit ball in this space is metrically convex.")
print("-" * 30)

# Define two points on the L-infinity unit sphere
u = np.array([1, 1])
v = np.array([-1, 1])

# Check that u and v are on the unit sphere
norm_u = l_inf_norm(u)
norm_v = l_inf_norm(v)

print(f"Point u = {list(u)}")
print(f"Norm of u: ||{list(u)}||_inf = max(|{u[0]}|, |{u[1]}|) = {norm_u}")
print(f"Point v = {list(v)}")
print(f"Norm of v: ||{list(v)}||_inf = max(|{v[0]}|, |{v[1]}|) = {norm_v}")
print("Both points are on the unit sphere.")
print("-" * 30)

# A point on the segment L(u,v) is given by (1-2t, 1) for t in [0,1].
# This is equivalent to (s, 1) for s in [-1,1].
# We verify that any such point has a norm of 1.
s_values = np.linspace(-1, 1, 5)
print("Verifying that the segment lies on the unit sphere:")
for s in s_values:
    point = np.array([s, 1])
    norm_point = l_inf_norm(point)
    print(f"For a point {np.round(point, 2).tolist()} on the segment, the norm is {norm_point}.")
print("The entire segment is on the unit sphere.")
print("-" * 30)

# Calculate the length of the segment, which is ||u - v||_inf
difference_vector = u - v
length = l_inf_norm(difference_vector)

print("Calculating the length of the line segment:")
final_equation = f"Length = ||u - v||_inf = ||{list(u)} - {list(v)}||_inf = ||{list(difference_vector)}||_inf = {length}"
print(final_equation)
print("\nConclusion: The largest possible length is 2.")
