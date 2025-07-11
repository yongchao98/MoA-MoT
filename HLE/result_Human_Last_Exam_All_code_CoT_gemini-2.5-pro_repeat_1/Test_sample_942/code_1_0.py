import numpy as np

# Step 1: Define the vector space and the functions.
# We use a discretized version of the space C[0,1].
# The variable 't' represents the domain [0,1].
t = np.linspace(0, 1, 1001)

# Define the functions c(t) and d(t) as numpy arrays.
c = np.ones_like(t)  # c(t) = 1
d = 1 - 2 * t      # d(t) = 1 - 2t

# The L-infinity norm is approximated by the maximum absolute value
# on the discretized domain.
def norm_inf(f):
    return np.max(np.abs(f))

# Step 2: Verify the properties of the functions c and d.
norm_c = norm_inf(c)
norm_d = norm_inf(d)
length = norm_inf(c - d)

print("Demonstration of a line segment of length 2 on the unit sphere.")
print("-" * 60)
print(f"The norm of function c(t)=1 is ||c|| = {norm_c}")
print(f"The norm of function d(t)=1-2t is ||d|| = {norm_d}")
print("Both functions are on the unit sphere.")
print("-" * 60)

# The "final equation" is the length of the segment.
print(f"The length of the segment is ||c - d||.")
# Output the numbers for the equation ||c-d|| = 2
print(f"Numerically calculated length: ||c - d|| = {length}")
print(f"Theoretical length: 2")
print("-" * 60)

# Step 3: Verify that the entire segment lies on the unit sphere.
# We check the norm for several points on the segment, indexed by s.
print("Verifying that segment points (1-s)c + sd have norm 1:")
all_norms_are_one = True
for s in np.linspace(0, 1, 11):
    segment_point = (1 - s) * c + s * d
    norm_s = norm_inf(segment_point)
    print(f"s = {s:.1f}, norm = {norm_s:.4f}")
    # Use a tolerance for floating point comparison
    if not np.isclose(norm_s, 1.0):
        all_norms_are_one = False

print("-" * 60)
if all_norms_are_one:
    print("Conclusion: The example is valid. A segment of length 2 exists.")
else:
    print("Conclusion: The numerical verification failed.")
