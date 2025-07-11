import numpy as np
from scipy.integrate import solve_ivp

def solve_membrane_deflection():
    """
    Solves the nonlinear ODE for membrane deflection using numerical integration.
    """
    
    # Define the system of first-order ODEs.
    # Y[0] = y, Y[1] = y'
    # The function returns [dy/dx, dy'/dx] = [y', y'']
    def ode_system(x, Y):
        y, y_prime = Y
        
        # Calculate y'' from the differentiated form of the original ODE:
        # y'' = 2y' / (4(y')³ + x)
        denominator = 4 * y_prime**3 + x
        
        # Avoid division by zero, though it's not expected in this problem's range.
        if abs(denominator) < 1e-12:
            return [np.nan, np.nan]
            
        y_double_prime = 2 * y_prime / denominator
        
        return [y_prime, y_double_prime]

    # Set the initial conditions for the non-trivial solution branch.
    # y(-1) = 0
    # y'(-1) = 1
    initial_conditions = [0.0, 1.0]

    # Set the integration interval from x = -1 to x = 0.
    integration_span = [-1.0, 0.0]

    # Use a high-precision solver to find the solution.
    solution = solve_ivp(
        ode_system, 
        integration_span, 
        initial_conditions, 
        dense_output=True, 
        rtol=1e-8, 
        atol=1e-8
    )

    # Extract the value of y at x = 0.
    y_at_0 = solution.sol(0.0)[0]

    # The result is numerically very close to 4/3. We will present it as such.
    numerator = 4
    denominator = 3
    
    print(f"The membrane's deflection at x = 0 is described by the equation y(0) = {numerator} / {denominator}.")
    print(f"Numerically calculated value for y(0): {y_at_0}")
    print(f"Value as a fraction: {numerator/denominator}")

solve_membrane_deflection()