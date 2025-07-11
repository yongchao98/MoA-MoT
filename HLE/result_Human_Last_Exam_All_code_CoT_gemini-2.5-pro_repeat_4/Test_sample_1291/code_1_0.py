# The problem asks to identify the types of points included in a specific subset of the Berkovich projective line.
# Based on the analysis, the construction maps pairs (z_0, z) from C_p^x * C_p to disks D(a, r)
# in the Berkovich line. The radius of these disks is r = 1/|z_0|_p.
# For any z_0 in C_p^x, its p-adic norm |z_0|_p must be in the value group of C_p, which is p^Q.
# This means the radius r is also in p^Q.
# By definition, points in the Berkovich line corresponding to disks with radii in the value group are Type 2 points.
# The construction cannot produce points of Type 1 (radius 0), Type 3 (radius not in p^Q), or Type 4 (limits of disks).
# Thus, the subset consists only of Type 2 points.

# The corresponding number for this type is 2.
print(2)