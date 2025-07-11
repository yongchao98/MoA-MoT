import numpy as np

# Based on the step-by-step derivation, the problem simplifies to finding
# the supremum of an angle, which is determined by the minimum possible
# cosine between the trajectory vector and the inner normal vector.

# 1. The relationship between the supremum angle M(theta) and the cosine is:
#    M(theta) = arccos(v_min)
#    where v_min is the minimum possible value of the dot product (cosine of the angle)
#    between the unit trajectory vector and the unit inner normal vector.

# 2. Our derivation showed that for any small theta > 0, it is possible to find
#    a trajectory that is perfectly anti-parallel to the inner normal vector.

# 3. When two vectors are anti-parallel, the cosine of the angle between them is -1.
#    Therefore, v_min = -1.

min_cosine_value = -1

# 4. We can now write the final equation for M(theta).
print("The supremum angle M(theta) is given by the equation:")
print(f"M(theta) = arccos({min_cosine_value})")

# 5. The value of M(theta) is constant for any small theta > 0.
M_theta = np.arccos(min_cosine_value)

# 6. The limit of a constant is the constant itself.
limit_M_theta = M_theta

print(f"The value of M(theta) is {M_theta}, which is pi.")
print(f"The limit of M(theta) as theta approaches 0 is therefore also pi.")

# The final numerical answer is pi.