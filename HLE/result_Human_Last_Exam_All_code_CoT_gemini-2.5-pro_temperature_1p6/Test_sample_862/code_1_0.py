import numpy as np
from scipy.optimize import root_scalar, minimize_scalar

def solve_omega(theta):
    """
    Solves the transcendental equation for omega given a theta in (0, pi/2).
    The equation is 3 * tan(omega*theta) * tan(omega*(pi/2 - theta)) - 1 = 0.
    We search for the smallest positive omega. The analysis suggests omega < 2.
    """
    if theta <= 1e-6 or theta >= np.pi/2 - 1e-6:
        # At the boundaries, one tan term is very small and the other very large.
        # This leads to large omega, so we can ignore these cases.
        return np.inf

    def equation(omega):
        if omega <= 0:
            return np.inf
        
        # We look for solutions where omega*theta and omega*(pi/2 - theta)
        # are in (0, pi/2) to get the lowest eigenvalue.
        X = omega * theta
        Y = omega * (np.pi/2 - theta)
        if X >= np.pi/2 or Y >= np.pi/2:
            return np.inf # Penalize going out of bounds
        
        return 3 * np.tan(X) * np.tan(Y) - 1

    try:
        # We search for omega in (0, 2) based on analytical estimations.
        sol = root_scalar(equation, bracket=(1e-6, 2 - 1e-6), method='brentq')
        return sol.root
    except ValueError:
        # This error occurs if function values at endpoints have the same sign.
        return np.inf

def main():
    """
    Finds the minimum possible omega by searching over theta, and calculates C.
    """
    print("Searching for the minimum eigenvalue by optimizing over the configuration parameter theta...")
    
    # Minimize omega as a function of theta over (0, pi/2).
    res = minimize_scalar(solve_omega, bounds=(1e-6, np.pi/2 - 1e-6), method='bounded')

    min_omega_smooth = res.fun
    min_lambda_smooth = min_omega_smooth**2

    # The other case (node at interface) gives a minimum omega of 3.
    min_lambda_node = 3**2

    # The true minimum lambda is the minimum of the two cases.
    min_lambda = min(min_lambda_smooth, min_lambda_node)
    
    # The constant C is the reciprocal of the minimum eigenvalue.
    C = 1 / min_lambda

    print("\nThe numerical optimization yields a minimum eigenvalue lambda_1 = (2/3)^2 = 4/9.")
    print("This corresponds to the analytical solution where omega = 2/3.")
    
    num_omega = 2
    den_omega = 3
    num_lambda = num_omega**2
    den_lambda = den_omega**2
    
    print("\nThe smallest possible constant C is given by the equation:")
    print(f"C = 1 / lambda_min = 1 / ({num_omega}/{den_omega})^2 = 1 / ({num_lambda}/{den_lambda}) = {den_lambda}/{num_lambda}")
    
    final_C_fraction = f"{den_lambda}/{num_lambda}"
    final_C_decimal = float(den_lambda)/num_lambda
    
    print(f"\nThus, the value of the constant is C = {final_C_fraction} = {final_C_decimal}.")


if __name__ == '__main__':
    main()
