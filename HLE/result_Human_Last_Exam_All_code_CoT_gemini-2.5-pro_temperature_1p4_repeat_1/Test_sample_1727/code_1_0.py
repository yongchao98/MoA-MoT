import numpy as np
from scipy.optimize import fsolve

def find_max_response_factor():
    """
    Calculates the stable amplitude factor chi by finding the maximum response
    of the system.
    """
    # Step 1: Find the value of H = kh that maximizes the response.
    # The condition for the maximum is tanh(H) = 1/H.
    # We define a function whose root is the solution.
    def func_to_solve(H):
        # We need to handle H=0 to avoid division by zero.
        if H == 0:
            return -np.inf
        return np.tanh(H) - 1/H

    # Use a numerical solver (fsolve) to find the root.
    # An initial guess of 1.0 is reasonable.
    H0 = fsolve(func_to_solve, 1.0)[0]

    # Step 2: Calculate chi using the derived simplified formula for the maximum response:
    # chi = -sqrt(H0^2 - 1) / H0^2
    # We will print out each component of this equation.

    print("The final equation for the stable amplitude factor is: chi = -sqrt(H0^2 - 1) / H0^2")
    print("Below are the values for each number in the equation:")
    
    # Calculate intermediate values
    H0_squared = H0**2
    H0_squared_minus_1 = H0_squared - 1
    sqrt_val = np.sqrt(H0_squared_minus_1)
    chi = -sqrt_val / H0_squared
    
    print(f"H0 (dimensionless wavenumber where response is maximum): {H0}")
    print(f"H0^2: {H0_squared}")
    print(f"H0^2 - 1: {H0_squared_minus_1}")
    print(f"sqrt(H0^2 - 1): {sqrt_val}")
    print(f"Final amplitude factor chi: {chi}")
    
    # Returning the final answer in the requested format
    return chi

if __name__ == '__main__':
    final_chi = find_max_response_factor()
    print(f"\n<<< {final_chi} >>>")
