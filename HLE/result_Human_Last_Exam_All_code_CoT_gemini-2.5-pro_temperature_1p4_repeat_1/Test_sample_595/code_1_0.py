import math

def max_crossings(dz):
    """
    Calculates the maximum number of grid lines crossed for a displacement dz.
    This is |dz| if dz is an integer, and floor(|dz|) + 1 otherwise.
    """
    if dz == int(dz):
        return abs(int(dz))
    else:
        return math.floor(abs(dz)) + 1

# The triangle has side lengths 18, 18, 18*sqrt(2).

# Case 1: Legs are aligned with the coordinate axes.
# The displacement vectors for the three sides are (18, 0), (0, 18), and (-18, 18).
s1_case1 = (18, 0)
s2_case1 = (0, 18)
s3_case1 = (-18, 18)

k1_s1 = max_crossings(s1_case1[0]) + max_crossings(s1_case1[1])
k1_s2 = max_crossings(s2_case1[0]) + max_crossings(s2_case1[1])
k1_s3 = max_crossings(s3_case1[0]) + max_crossings(s3_case1[1])
k1_total = k1_s1 + k1_s2 + k1_s3

# Case 2: Hypotenuse is aligned with a coordinate axis.
# The legs are at 45 degrees to the axes. The displacement vectors are:
# Leg 1: (18/sqrt(2), 18/sqrt(2)) = (9*sqrt(2), 9*sqrt(2))
# Leg 2: (9*sqrt(2), -9*sqrt(2))
# Hypotenuse: (-18*sqrt(2), 0)
leg_disp = 9 * math.sqrt(2)
hyp_disp = 18 * math.sqrt(2)

s1_case2 = (leg_disp, leg_disp)
s2_case2 = (leg_disp, -leg_disp)
s3_case2 = (-hyp_disp, 0)

k2_s1 = max_crossings(s1_case2[0]) + max_crossings(s1_case2[1])
k2_s2 = max_crossings(s2_case2[0]) + max_crossings(s2_case2[1])
k2_s3 = max_crossings(s3_case2[0]) + max_crossings(s3_case2[1])
k2_total = k2_s1 + k2_s2 + k2_s3

# The largest k is the answer.
# We will print the calculation for the optimal case.
if k1_total > k2_total:
    print("Optimal orientation: Legs aligned with axes.")
    print(f"Number of squares crossed by side 1: {k1_s1}")
    print(f"Number of squares crossed by side 2: {k1_s2}")
    print(f"Number of squares crossed by side 3: {k1_s3}")
    print(f"Total squares crossed: {k1_s1} + {k1_s2} + {k1_s3} = {k1_total}")
    final_k = k1_total
else:
    print("Optimal orientation: Hypotenuse aligned with an axis.")
    print(f"Number of squares crossed by side 1 (leg): {k2_s1}")
    print(f"Number of squares crossed by side 2 (leg): {k2_s2}")
    print(f"Number of squares crossed by side 3 (hypotenuse): {k2_s3}")
    print(f"Total squares crossed: {k2_s1} + {k2_s2} + {k2_s3} = {k2_total}")
    final_k = k2_total
