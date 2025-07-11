import sympy as sp

# Step 1: Define the symbols for our coordinates
r, theta, z = sp.symbols('r theta z')
pi = sp.pi

# Step 2: Define the function and the integrand in cylindrical coordinates
# The function is f = z^2 * (x^2 + y^2) = z^2 * r^2
# The volume element dV is r*dr*d(theta)*dz
# The full term to integrate is (z^2 * r^2) * r
integrand = z**2 * r**3

# Step 3: Define the limits of integration
# z goes from 0 to 2
z_limits = (z, 0, 2)
# theta goes from 0 to 2*pi
theta_limits = (theta, 0, 2*pi)
# r goes from 0 to the edge of the cone, which is 3*(1 - z/2)
r_limits = (r, 0, 3 * (1 - z/2))

# Step 4: Compute the integral
# We integrate from the inside out: first r, then theta, then z.
integral_result = sp.integrate(integrand, r_limits, theta_limits, z_limits)

# Step 5: Display the results
print("Calculating the integral of f(x,y,z) = z^2*(x^2+y^2) over the specified cone.")
print(f"The integral is set up in cylindrical coordinates as: ∫(z=0 to 2) ∫(θ=0 to 2π) ∫(r=0 to 3(1-z/2)) {integrand} dr dθ dz")

# The result from SymPy is a symbolic expression of the form (numerator * pi / denominator)
# We can extract these parts to display them clearly.
if isinstance(integral_result, sp.Mul) and pi in integral_result.args:
    coefficient = integral_result.as_coeff_Mul(pi)[0]
    numerator = sp.numer(coefficient)
    denominator = sp.denom(coefficient)
    print("\n--- Result ---")
    print(f"The exact value of the integral is ({numerator} * π) / {denominator}")
    
    # As requested, output the numbers in the final equation
    print(f"\nFinal Equation: Integral = {numerator} * π / {denominator}")
    
    numerical_value = float(integral_result)
    print(f"The numerical value is approximately: {numerical_value:.4f}")
else:
    # Fallback for unexpected result format
    print(f"The value of the integral is: {integral_result}")
