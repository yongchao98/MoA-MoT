import numpy as np

def linf_norm(vector):
    """Calculates the L-infinity norm of a vector."""
    return np.max(np.abs(vector))

# Define two points u and v on the L_infinity unit sphere in R^2.
# These points are the endpoints of one of the sides of the square unit sphere.
u = np.array([1.0, -1.0])
v = np.array([1.0, 1.0])

print(f"Consider the vector space R^2 with the L_infinity norm (max(|x|, |y|)).")
print(f"Let u = {u} and v = {v}.")
print(f"The L_infinity norm of u is ||u||_inf = {linf_norm(u)}.")
print(f"The L_infinity norm of v is ||v||_inf = {linf_norm(v)}.")
print("Both points lie on the unit sphere.")
print("\nLet's check if the line segment between u and v lies on the unit sphere.")
print("A point x on the segment is given by x = (1-t)*u + t*v for t in [0,1].")

is_on_sphere = True
# Check 11 points on the line segment
for t in np.linspace(0, 1, 11):
    x = (1 - t) * u + t * v
    norm_x = linf_norm(x)
    print(f"For t = {t:.1f}, x = {x}, and ||x||_inf = {norm_x:.2f}")
    # Using a tolerance for floating point comparison
    if not np.isclose(norm_x, 1.0):
        is_on_sphere = False

if is_on_sphere:
    print("\nConclusion: The entire line segment lies on the unit sphere.")
else:
    print("\nError: The line segment does not lie on the unit sphere.")

# Calculate the length of the segment, which is the norm of the difference vector
length = linf_norm(u - v)

print("\nThe length of the line segment is the norm of the difference vector u-v.")
print(f"The final equation for the length is:")
print(f"||{u} - {v}||_inf = {length}")
