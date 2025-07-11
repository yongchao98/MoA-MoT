import sympy

# Set up the symbols
pi = sympy.pi
r, theta, z = sympy.symbols('r theta z')

# Explain the problem and the coordinate system change
print("To find the integral of f(x,y,z) = z^2(x^2+y^2) over the given cone, we use cylindrical coordinates.")
print("The transformation is x = r*cos(theta), y = r*sin(theta), z = z.")
print("The function becomes f(r,theta,z) = z^2 * r^2.")
print("The volume element is dV = r dr dtheta dz.")
print("So, the expression to integrate is (z^2 * r^2) * r = z^2 * r^3.")
print("")

# Define the integrand and limits
integrand = z**2 * r**3
# The equation for the cone surface is r = (R/H)*(H-z) = (3/2)*(2-z)
r_upper_limit = sympy.Rational(3, 2) * (2 - z)

# Explain the integral setup
print("The integral is set up as:")
print(f"I = Integral from 0 to 2*pi [ Integral from 0 to 2 [ Integral from 0 to {r_upper_limit} [ {integrand} ] dr ] dz ] dtheta")
print("")

# Perform and explain each step of the integration
# Step 1: Integrate with respect to r
print("Step 1: Integrating with respect to r...")
integral_r = sympy.integrate(integrand, (r, 0, r_upper_limit))
print(f"The integral of {integrand} from r=0 to r={r_upper_limit} is: {integral_r}")
print("")

# Step 2: Integrate with respect to z
print("Step 2: Integrating the result with respect to z...")
integral_z = sympy.integrate(integral_r, (z, 0, 2))
print(f"The integral of ({integral_r}) from z=0 to z=2 is: {integral_z}")
print("")

# Step 3: Integrate with respect to theta
print("Step 3: Integrating the result with respect to theta...")
final_integral = sympy.integrate(integral_z, (theta, 0, 2 * pi))
print(f"The integral of {integral_z} from theta=0 to theta=2*pi is: {final_integral}")
print("")

# Display the final answer clearly
print("The final equation for the integral is:")
coeff = final_integral.coeff(pi)
num, den = coeff.as_numer_denom()
print(f"I = ({num} * pi) / {den}")

# The final numeric value can be evaluated as:
# print(f"I = {final_integral.evalf()}")
<<<108*pi/35>>>