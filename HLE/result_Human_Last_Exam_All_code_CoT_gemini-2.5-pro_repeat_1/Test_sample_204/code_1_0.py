import sympy

# Define symbols
R = sympy.Symbol('R', real=True, positive=True)
pi = sympy.pi

# The Hopf charge for this axially symmetric field configuration is given by
# the formula Q = -k * N.
# k is the winding number of the azimuthal angle 'f'.
# N is the degree (topological charge) of the profile map in the poloidal plane.

# 1. Determine the azimuthal winding number, k.
# The field's azimuthal angle is f = atan2(y, x). This is the standard
# cylindrical angle phi. As one circles the z-axis, f sweeps from 0 to 2*pi.
# Therefore, the winding number k is 1.
k = 1

# 2. Determine the degree of the profile map, N.
# The profile is determined by the polar angle G, which depends on
# r2 = sqrt((x*x+y*y-0.5)*(x*x+y*y-0.5)+z*z).
# This structure corresponds to a 2D skyrmion in a plane defined by coordinates
# u = x*x+y*y-0.5 and v=z.
# The radial coordinate in this plane is R = sqrt(u^2+v^2), which is exactly r2.
# The degree N is calculated by the integral:
# N = (1 / (4*pi)) * integral( sin(G) dG wedge dTheta ) over the (u,v) plane.
# This simplifies to N = (1/2) * integral_from_0_to_inf( sin(G(R)) * dG/dR ) dR.

# Define G as a function of R
G_of_R = pi * sympy.exp(-10 * R)

# Calculate the derivative dG/dR
dG_dR = sympy.diff(G_of_R, R)

# The integrand is sin(G(R)) * dG/dR
integrand = sympy.sin(G_of_R) * dG_dR

# Perform the symbolic integration from R=0 to R=infinity
# This corresponds to integrating over the entire (u,v) plane.
# The limits for G are: R=0 -> G=pi; R=inf -> G=0.
integral_val = sympy.integrate(integrand, (R, 0, sympy.oo))

# The degree N is half of the integral value.
N = integral_val / 2

# 3. Calculate the Hopf charge Q.
Q = -k * N

# Print the results of each step and the final equation.
print(f"The azimuthal winding number is k = {k}")
print(f"The degree of the profile map is N = {N}")
print(f"The Hopf charge is given by the equation Q = -k * N")
print(f"Substituting the values, we get:")
print(f"Q = -({k}) * ({N})")
final_Q = -k * N
print(f"Q = {final_Q}")

<<<1>>>