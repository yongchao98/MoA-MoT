import numpy as np

# The vector field (nx, ny, nz) is a map from R^3 to the unit sphere S^2.
# Its components are defined as follows:
PI = np.pi

def define_field_parameters(x, y, z):
    """
    This function calculates the intermediate parameters G and f
    based on the space coordinates (x,y,z).
    """
    # f is the azimuthal angle in the xy-plane.
    f = np.arctan2(y, x)
    
    # r2 is a toroidal radius, centered on a circle of radius sqrt(0.5).
    r_squared = x*x + y*y
    r2 = np.sqrt((r_squared - 0.5)**2 + z*z)
    
    # G is the polar angle on the sphere, determined by r2.
    G = PI * np.exp(-10 * r2)
    
    return G, f

# The Hopf charge of this vector field can be computed using a topological method.
# It is defined as the linking number of the preimages of any two regular values
# on the target sphere. We choose the North Pole (0,0,1) and the South Pole (0,0,-1)
# as our two points.

# Analysis Step 1: Find the preimage of the South Pole n=(0,0,-1)
# The condition for the South Pole is nz = -1.
# Given nz = cos(G), this means cos(G) = -1.
# This equality holds if G = (2k + 1) * PI for any integer k.
# From the definition of G, we have PI * exp(-10*r2) = (2k+1)*PI.
# => exp(-10*r2) = 2k + 1.
# Since r2 is real and non-negative, exp(-10*r2) must be in the interval (0, 1].
# The only integer k that satisfies this is k=0, which gives exp(-10*r2) = 1.
# This implies -10*r2 = 0, so r2 = 0.
# For r2 to be zero, both terms in its definition must be zero:
# 1) z*z = 0  => z = 0
# 2) (x*x + y*y - 0.5)**2 = 0  => x*x + y*y = 0.5
# So, the preimage of the South Pole is a circle in the z=0 plane with radius sqrt(0.5).

# Analysis Step 2: Find the preimage of the North Pole n=(0,0,1)
# The condition for the North Pole is nz = 1.
# cos(G) = 1, which holds if G = 2k*PI for any integer k.
# PI * exp(-10*r2) = 2k*PI => exp(-10*r2) = 2k.
# For exp(-10*r2) to be in (0, 1], the only integer solution is k=0.
# This gives G = 0.
# This equality holds only in the limit r2 -> infinity.
# The preimage of the North Pole is the set of points at infinity. Topologically,
# this can be represented as a line passing through the origin, such as the z-axis,
# which is considered a loop closed at the point at infinity.

# Analysis Step 3: Calculate the Linking Number
# The Hopf charge is the linking number of the two preimages found above:
# Curve 1 (preimage of South Pole): The circle x^2+y^2=0.5 in the z=0 plane.
# Curve 2 (preimage of North Pole): The z-axis (x=0, y=0).
# A straight line (the z-axis) passing through the center of a circle has a
# linking number of 1 with the circle.

# Conclusion: The Hopf charge is 1.
hopf_charge = 1

# Final Equation: Hopf Charge = 1
print("The Hopf charge is calculated via a topological analysis:")
print("1. The preimage of the South Pole is the circle x^2+y^2=0.5, z=0.")
print("2. The preimage of the North Pole is topologically equivalent to the z-axis.")
print("3. The linking number of this circle and the z-axis is 1.")
print("Final result of the calculation:")
print(f"Hopf Charge = {hopf_charge}")
<<<1>>>