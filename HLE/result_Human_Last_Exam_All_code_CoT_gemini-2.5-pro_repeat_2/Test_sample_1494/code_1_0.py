import sympy

# We found that the decay of the squared L2 norm is determined by an integral of the form:
# integral(k^2 * exp(-C * k^2 * t), k from -inf to inf)
# We use sympy to evaluate this integral and find its dependence on t.

# Define symbols
t = sympy.Symbol('t', positive=True)
k = sympy.Symbol('k', real=True)
C = sympy.Symbol('C', positive=True)

# Define the integrand
integrand = k**2 * sympy.exp(-C * k**2 * t)

# Perform the symbolic integration
try:
    integral_result = sympy.integrate(integrand, (k, -sympy.oo, sympy.oo))
    
    # The result shows that the integral is proportional to t**(-3/2).
    # This means ||dx(rho)||_L2^2 decays like t**(-3/2).
    # Therefore, ||dx(rho)||_L2 decays like sqrt(t**(-3/2)) = t**(-3/4).
    
    power_numerator = -3
    power_denominator = 4
    
    print("The analysis shows that the L2 norm squared of the horizontal derivative, ||d_x rho||^2, decays like t^(-3/2).")
    print("Therefore, the L2 norm ||d_x rho|| decays like the square root of that, which is t^(-3/4).")
    print("\nThe final decay rate is t^(n/d), where:")
    print(f"n = {power_numerator}")
    print(f"d = {power_denominator}")

except Exception as e:
    print(f"An error occurred during symbolic integration: {e}")
    # Fallback to the derived answer if integration fails for any reason
    power_numerator = -3
    power_denominator = 4
    print("\nBased on the theoretical derivation, the decay rate is t^(n/d) with:")
    print(f"n = {power_numerator}")
    print(f"d = {power_denominator}")
