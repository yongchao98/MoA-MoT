import sympy as sp

# Step 1: Define symbolic parameters for the vehicle model
# a: distance from CG to front axle
# b: distance from CG to rear axle
# c_f: cornering stiffness of front axle
# c_r: cornering stiffness of rear axle
# m: vehicle mass
# I: vehicle moment of inertia
# v: vehicle forward speed
a, b, c_f, c_r, m, I, v = sp.symbols('a b c_f c_r m I v', positive=True)

print("Deriving the critical speed for an oversteering vehicle.")
print("Step 1: Define the state matrix A for the linear single-track model.")
print("The state vector is [sideslip_angle, yaw_rate]^T.\n")

# Step 2: Define the elements of the state matrix A
# The matrix is derived from the vehicle's lateral and yaw dynamics equations.
A11 = -(c_f + c_r) / (m * v)
A12 = (b * c_r - a * c_f) / (m * v**2) - 1
A21 = (b * c_r - a * c_f) / I
A22 = -(a**2 * c_f + b**2 * c_r) / (I * v)

A = sp.Matrix([[A11, A12], [A21, A22]])

print("State Matrix A:")
sp.pprint(A)
print("\n" + "="*50 + "\n")

print("Step 2: Analyze stability using the characteristic equation det(s*I - A) = 0.")
print("The equation is s^2 + k1*s + k0 = 0.")
print("Stability is lost when the constant term k0 becomes zero.\n")

# Step 3: Calculate the constant term k0, which is equal to the determinant of A.
k0 = sp.simplify(sp.det(A))

print("The constant term k0 of the characteristic polynomial is:")
sp.pprint(k0)
print("\n" + "="*50 + "\n")

print("Step 4: Solve for the speed v where k0 = 0 to find the critical speed.")
print("For an oversteering vehicle, (a*c_f - b*c_r) is positive.")
print("Setting k0 = 0 and solving for v^2...\n")

# Step 5: Solve the equation k0 = 0 for v^2.
# The solver will find the value of v^2 that makes the numerator of k0 zero.
v_crit_squared_sol = sp.solve(sp.Eq(k0, 0), v**2)

# The solution is a list containing one element.
v_crit_squared = v_crit_squared_sol[0]

# Extract the numerator and denominator of the final expression for clarity.
numerator = sp.numer(v_crit_squared)
denominator = sp.denom(v_crit_squared)

print("The expression for the square of the critical speed (v_crit^2) is:")
print("\nv_crit^2 = (", numerator, ") / (", denominator, ")")
print("-" * 50)

# Final step: Take the square root to find the critical speed v_crit.
v_crit = sp.sqrt(v_crit_squared)

print("\nThe final formula for the critical speed, v_crit, is:")
sp.pprint(v_crit)

# The final result is converted to a string for the required output format.
final_answer = str(v_crit)
print(f"\n<<<{final_answer}>>>")