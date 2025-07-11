import sympy

# 1. Define symbolic variables for the vehicle parameters.
# The parameters are assumed to be positive real numbers.
# a: distance from CG to front axle
# b: distance from CG to rear axle
# m: vehicle mass
# I: vehicle moment of inertia around the vertical axis
# c_f: front axle cornering stiffness
# c_r: rear axle cornering stiffness
# v: forward speed
a, b, m, I, c_f, c_r, v = sympy.symbols('a b m I c_f c_r v', positive=True)

# 2. Define the state matrix A for the linear single-track model.
# The state vector is x = [v_y, r]^T, where v_y is lateral velocity and r is yaw rate.
# The system is described by ẋ = Ax. We analyze stability by examining A.

A11 = -(c_f + c_r) / (m * v)
A12 = -v - (a * c_f - b * c_r) / (m * v)
A21 = (-a * c_f + b * c_r) / (I * v)
A22 = -(a**2 * c_f + b**2 * c_r) / (I * v)

A = sympy.Matrix([[A11, A12], [A21, A22]])

# 3. Calculate the determinant of the state matrix A.
# The system becomes unstable when the determinant equals zero.
det_A = sympy.det(A)
det_A_simplified = sympy.simplify(det_A)

# 4. Create the equation det(A) = 0 to find the critical speed.
# The critical speed, v_crit, is the speed 'v' that satisfies this equation.
# Note: For a real solution for v_crit to exist, the vehicle must be
# oversteering, meaning a*c_f > b*c_r.
stability_equation = sympy.Eq(det_A_simplified, 0)

# 5. Solve the equation for the critical speed v.
# Sympy's solve function will find the expression for 'v'.
solutions = sympy.solve(stability_equation, v)

# The result is a list of solutions. We are interested in the positive one.
v_crit = solutions[0]

# 6. Print the derived formula for the critical speed.
# This code block outputs each component of the final equation symbolically.
print("The derived formula for the critical speed (v_crit) is:")
equation_to_print = sympy.Eq(sympy.Symbol('v_crit'), v_crit)
sympy.pretty_print(equation_to_print)
