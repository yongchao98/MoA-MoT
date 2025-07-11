import sympy as sp

# Step 1 & 2: Define the geometry and the goal in symbolic terms
theta = sp.Symbol('theta', real=True, positive=True)

# Vertices and side A
# O = (0, 0)
# P = (5, 0)
Q_x, Q_y = 5 * sp.cos(theta), 5 * sp.sin(theta) # Point Q

# Unit circle arc endpoints
x1_x, x1_y = 1, 0                     # Point x1 (at angle 0)
x2_x, x2_y = sp.cos(theta), sp.sin(theta) # Point x2 (at angle theta)

# Step 3: Find the range of trajectory directions for small theta
# We analyze the directions of the lines connecting the endpoints of the domains.
# Trajectory 1: from x1 to Q
v1_x = Q_x - x1_x
v1_y = Q_y - x1_y
# Angle of trajectory 1, for small theta. Slope is v1_y / v1_x.
# tan(psi1) approx (5*theta) / (5-1) = 5*theta/4. So psi1 approx 5*theta/4
psi1 = sp.atan2(v1_y, v1_x)

# Trajectory 2: from x2 to P
v2_x = 5 - x2_x
v2_y = 0 - x2_y
# Angle of trajectory 2, for small theta. Slope is v2_y / v2_x
# tan(psi2) approx -theta / (5-1) = -theta/4. So psi2 approx -theta/4
psi2 = sp.atan2(v2_y, v2_x)

# For small theta, the range of trajectory angles psi is [psi2, psi1]
# We will use the small angle approximations: psi_min = -theta/4, psi_max = 5*theta/4

# Step 4: Find the direction of the inner normal to side A
# Side A connects P and Q. Vector along A is vec_A = Q - P
vec_A_x = Q_x - 5
vec_A_y = Q_y - 0
# A normal vector is perpendicular to vec_A.
# Outer normal points away from the origin. For small theta, Q is above P, so
# outer normal points right and up. vec_A is approx (0, 5*theta).
# So normal is approx (1,0) direction. More accurately, angle is theta/2.
# Outer normal direction vector:
n_out_x, n_out_y = sp.cos(theta/2), sp.sin(theta/2)
# Inner normal is in the opposite direction (angle + pi)
n_in_x, n_in_y = -sp.cos(theta/2), -sp.sin(theta/2)
# Inner normal vector direction
u_norm = sp.Matrix([n_in_x, n_in_y])

# Step 5: Find the supremum angle M(theta)
# The angle alpha is given by cos(alpha) = |u_traj . u_norm|
# where u_traj is the direction vector of the trajectory (cos(psi), sin(psi))
# cos(alpha) = |cos(psi)*n_in_x + sin(psi)*n_in_y|
# = |cos(psi)*(-cos(theta/2)) + sin(psi)*(-sin(theta/2))|
# = |-cos(psi - theta/2)| = |cos(psi - theta/2)|
#
# To maximize alpha, we must minimize cos(alpha)
# We need to find min |cos(psi - theta/2)| for psi in [-theta/4, 5*theta/4]
# Let z = psi - theta/2.
# The range for z is [-theta/4 - theta/2, 5*theta/4 - theta/2] = [-3*theta/4, 3*theta/4]
# For small theta, z is small, so cos(z) is positive.
# We need to minimize cos(z) on z in [-3*theta/4, 3*theta/4].
# cos(z) is an even function and decreasing for z in [0, pi].
# The minimum occurs at the endpoints of the interval, z = +/- 3*theta/4.
# min(cos(z)) = cos(3*theta/4).
# So, cos(M(theta)) = cos(3*theta/4).
# Since M(theta) is the angle of incidence, it's in [0, pi/2].
# For small theta, 3*theta/4 is also in this range.
# Therefore, M(theta) = 3*theta/4.

# Step 6: Find the limit of M(theta) as theta -> 0
# M(theta) is a function of theta. Let's define it based on our finding.
M_theta = (sp.S(3)/sp.S(4)) * theta

# Now compute the limit
limit_M = sp.limit(M_theta, theta, 0)

print("For small theta, M(theta) is approximately 3*theta/4.")
# The question asks for the limit value of M(theta)
print("The limit of M(theta) as theta goes to 0 is:")
# The sp.limit will return a sympy object, we cast it to an int for cleaner output.
print(int(limit_M))
