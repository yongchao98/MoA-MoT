import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit

def power_law(t, a, b):
    """Defines a power law relationship."""
    return a * t**b

def calculate_norm_squared(t, C=1.0, initial_amplitude=1.0):
    """
    Calculates the integral representing the squared L2 norm of the horizontal derivative.
    The integral is over k from 0 to infinity, multiplied by 2 to account for negative k.
    We assume |rho_hat_0(k)|^2 is constant (initial_amplitude^2) for small k.
    """
    # Integrand is k^2 * |rho_hat_0|^2 * exp(-2*t*sigma_k)
    # For small k, sigma_k ~ C * k^2
    integrand = lambda k: (k**2) * (initial_amplitude**2) * np.exp(-2 * t * C * k**2)
    # The analytical value of the integral from 0 to infinity is (sqrt(pi)/4) * (2*t*C)^(-3/2)
    # We compute it numerically to demonstrate.
    result, _ = quad(integrand, 0, np.inf)
    return 2 * result # Factor of 2 for symmetry k and -k

def main():
    """
    Main function to simulate decay and fit a power law.
    """
    # Time points for calculation
    t_values = np.logspace(0, 4, 20)
    norm_sq_values = []

    print("Calculating norm squared at different time points...")
    for t in t_values:
        norm_sq = calculate_norm_squared(t, C=0.1)
        norm_sq_values.append(norm_sq)
        # print(f"t = {t:.2e}, ||dx rho||^2 = {norm_sq:.2e}")

    norm_sq_values = np.array(norm_sq_values)

    # Fit a power law to the calculated values
    try:
        # We fit log(y) = log(a) + b*log(x) which is linear
        log_t = np.log(t_values)
        log_norm_sq = np.log(norm_sq_values)
        params, covariance = np.polyfit(log_t, log_norm_sq, 1, cov=True)
        slope = params[0]
        intercept = params[1]
        
        print("\nFitting log(||dx rho||^2) vs log(t) to a line 'y = slope*x + intercept'")
        print(f"Fit results:")
        print(f"  Slope (Decay Exponent for norm squared): {slope:.4f}")
        print(f"  Intercept: {intercept:.4f}")
        
        decay_exponent = slope / 2
        print(f"\nThe decay rate for the norm ||dx rho(t)||_L2 is t^{slope:.4f}/2 = t^{decay_exponent:.4f}.")
        print("This is very close to the theoretical value of t^(-3/4) = t^(-0.75).")

    except Exception as e:
        print(f"Could not perform fit: {e}")

    # Theoretical integral value for k^2*exp(-a*k^2) is sqrt(pi)/(4*a^(3/2))
    # where a = 2*t*C. So result should be proportional to t^(-3/2).
    # This implies the slope of log-log plot is -1.5
    print("\nTheoretical decay exponent for the squared norm is -1.5.")
    
    # Let's print the theoretical result for comparison
    C_const = 0.1
    initial_amplitude = 1.0
    a = 2 * t_values * C_const
    theoretical_integral_analytic = (np.sqrt(np.pi) / (4 * a**(3/2)))
    
    # We should have norm_sq_values close to 2 * initial_amplitude**2 * theoretical_integral_analytic
    # which indeed scales as t^(-3/2).
    
if __name__ == '__main__':
    main()
