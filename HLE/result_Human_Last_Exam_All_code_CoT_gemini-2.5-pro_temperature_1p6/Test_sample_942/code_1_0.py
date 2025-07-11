import numpy as np

def linf_norm(v):
    """Calculates the L-infinity norm of a vector."""
    return np.max(np.abs(v))

# Part 1: Define two points u and v in R^2
# In a space where the unit ball is metrically convex (e.g., R^n with L_infinity norm),
# we can find a line segment on the unit sphere with length 2.
u = np.array([1.0, -1.0])
v = np.array([1.0, 1.0])

# Part 2: Verify that u and v are on the unit sphere
u_norm = linf_norm(u)
v_norm = linf_norm(v)

print(f"Vector u = {u}, Norm ||u||_inf = {u_norm}")
print(f"Vector v = {v}, Norm ||v||_inf = {v_norm}")
if u_norm == 1 and v_norm == 1:
    print("Both u and v are on the unit sphere.")
else:
    print("Error: u or v are not on the unit sphere.")

# Part 3: Verify that the line segment between u and v is on the unit sphere
# We check a few points on the segment L(u,v) = (1-t)u + tv for t in [0,1]
segment_on_sphere = True
print("\nChecking points on the segment L(u,v):")
for t in np.linspace(0, 1, 11):
    point = (1 - t) * u + t * v
    norm_point = linf_norm(point)
    # Using np.isclose for robust float comparison
    if not np.isclose(norm_point, 1.0):
        segment_on_sphere = False
    print(f"t={t:.1f}, point={np.round(point, 2)}, norm={norm_point:.2f}")

if segment_on_sphere:
    print("\nAll tested points are on the unit sphere. The segment lies on the unit sphere.")
else:
    print("\nThe segment does not lie on the unit sphere.")

# Part 4: Calculate the length of the segment, which is ||u-v||
diff_vector = u - v
length = linf_norm(diff_vector)

print("\nThe length of the segment is the norm of the difference vector.")
print(f"u - v = ({u[0]},{u[1]}) - ({v[0]},{v[1]}) = ({diff_vector[0]},{diff_vector[1]})")
print(f"||u - v||_inf = ||({diff_vector[0]}, {diff_vector[1]})||_inf = {length}")

print("\nAs shown by derivation, the maximum possible length is 2, and this example achieves it.")
print("The final equation for the length calculation is:")
# Output the final equation as requested
print(f"||({u[0]}-{v[0]}, {u[1]}-{v[1]})|| = ||({diff_vector[0]}, {diff_vector[1]})|| = {length}")
