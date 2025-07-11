import numpy as np

# The feasible set F is the line segment { (2, 0, x3) | -2 <= x3 <= -1 }.
# The point x* is (2, 0, -1), an endpoint of the segment.
# The tangent cone T_F(x*) consists of vectors d = (d1, d2, d3) pointing from x* into F.
# This gives T_F(x*) = {d in R^3 | d1=0, d2=0, d3 <= 0}.
# The normal cone T_F°(x*) is the polar of the tangent cone.
# A vector s = (s1, s2, s3) is in the normal cone if s^T * d <= 0 for all d in T_F(x*).
# s^T * d = s1*d1 + s2*d2 + s3*d3 = s3*d3.
# The condition is s3*d3 <= 0 for all d3 <= 0. This implies s3 >= 0.
# There are no constraints on s1 and s2.

print("An explicit representation of the normal cone T_F°(x*) is the set of vectors s = (s_1, s_2, s_3) in R^3 that satisfy the following conditions:")
print("s_1: can be any real number")
print("s_2: can be any real number")
print("s_3: must be greater than or equal to 0")
print("\nIn set notation, the normal cone is:")
print("T_F°(x*) = { s = (s_1, s_2, s_3) ∈ R^3 | s_3 >= 0 }")

print("\nThe final equation that defines the cone is:")
s_index = 3
inequality_symbol = ">="
value = 0
print(f"s_{s_index} {inequality_symbol} {value}")
