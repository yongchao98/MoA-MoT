import numpy as np
from scipy.integrate import quad

def solve_ode_approximation():
    """
    This function calculates the parameters for the approximate analytical solution
    to the given ODE in the large x regime (i.e., near its singularity).
    """
    # Step 1: Define initial conditions for u = y'.
    # From the problem statement, y'(0) = 3.00 and y''(0) = 2.00.
    u0 = 3.00
    u_prime_0 = 2.00

    # Step 2: Use the simplified ODE u'' = u^4 to find the first integral.
    # 0.5 * (u')^2 = 0.2 * u^5 + K
    # Calculate the integration constant K using initial conditions u(0) and u'(0).
    K = 0.5 * u_prime_0**2 - 0.2 * u0**5

    # Step 3: From the integrated equation, u' = sqrt(0.4*u^5 + 2*K).
    # We can now set up integrals to find the singularity time x_s and the value y(x_s) = C0.
    
    # Integrand for calculating x_s, which is the integral of dt = du/u'
    def integrand_xs(u):
        value_under_sqrt = 0.4 * u**5 + 2 * K
        if value_under_sqrt < 0:
            return 0 # The integration range should not result in this.
        return 1.0 / np.sqrt(value_under_sqrt)

    # Integrand for calculating C0 = y(x_s). C0 is the integral of u(t)dt = u*du/u'
    def integrand_C0(u):
        value_under_sqrt = 0.4 * u**5 + 2 * K
        if value_under_sqrt < 0:
            return 0
        return u / np.sqrt(value_under_sqrt)

    # Step 4: Calculate x_s and C0 using numerical integration (quad).
    # The integration is from the initial value u0 to infinity.
    xs, _ = quad(integrand_xs, u0, np.inf)
    
    # The constant C0 in the final solution y(x) = C0 - C1*(xs-x)^(1/3) is y(xs).
    # y(xs) = y(0) + integral from 0 to xs of y'(t) dt = 0 + integral of u(t)dt.
    C0, _ = quad(integrand_C0, u0, np.inf)

    # Step 5: Calculate the coefficient C1.
    # From the dominant balance analysis, C1 = 3 * (10/9)^(1/3).
    C1 = 3 * (10/9)**(1/3)

    # Step 6: Print the final analytical expression with coefficients rounded to two decimal places.
    print("The approximate analytical expression for y(x) in the large x regime is:")
    # The final equation includes the numbers for C0, C1, and xs, all formatted to two decimal places.
    print(f"y(x) = {C0:.2f} - {C1:.2f} * ({xs:.2f} - x)^(1/3)")

solve_ode_approximation()