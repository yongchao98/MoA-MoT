import numpy as np
from scipy.integrate import quad

def main():
    """
    This script estimates the coefficient C for the scaling law of the fluctuation R.
    The scaling law is R(epsilon) = C * epsilon^0.5.
    """
    
    # The Green's function G(x,s) for the operator L[u] = u_xx - u_x on [0,1]
    # with boundary conditions u(0)=u(1)=0 is derived using standard methods.
    # It solves L[G] = delta(x-s).
    def G(x, s):
        e = np.e
        # The Wronskian for the fundamental solutions u1=e^x-1 and u2=e^x-e
        w_s = (e - 1) * np.exp(s)
        
        # Fundamental solutions satisfying one boundary condition each
        u1_x = np.exp(x) - 1
        u2_x = np.exp(x) - e
        u1_s = np.exp(s) - 1
        u2_s = np.exp(s) - e
        
        if x < s:
            return (u1_x * u2_s) / w_s
        else: # x >= s
            return (u1_s * u2_x) / w_s

    # Function to compute the variance of G(x,s) for a fixed x, where s ~ Uniform(0,1)
    def get_variance_G(x):
        # E[G] = integral of G(x,s) ds from 0 to 1
        E_G, _ = quad(lambda s: G(x, s), 0, 1, limit=100)
        # E[G^2] = integral of G(x,s)^2 ds from 0 to 1
        E_G_sq, _ = quad(lambda s: G(x, s)**2, 0, 1, limit=100)
        
        return E_G_sq - E_G**2

    # Find the maximum of the variance over x in [0,1].
    # We sample x in (0,1) to avoid potential issues at the boundaries.
    X_values = np.linspace(0.01, 0.99, 300)
    variances = np.array([get_variance_G(x) for x in X_values])
    
    # The constant C is the square root of the maximum variance.
    max_variance = np.max(variances)
    C = np.sqrt(max_variance)

    print("The estimated maximum magnitude of fluctuations is R(epsilon) = C * epsilon^p.")
    print("\nThe values for the parameters in this equation are:")
    # Print each number in the final equation as requested.
    print(f"The coefficient C is: {C:.6f}")
    print(f"The scaling power p is: {0.5}")

if __name__ == '__main__':
    main()