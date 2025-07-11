import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_expression():
    """
    Solves the ODE numerically to find the constants for the asymptotic expression.
    """
    # The ODE is approximated by replacing the oscillating tan term with its average value 0.5.
    # This simplifies numerical integration and is justified because the (y')^4 term
    # dominates near the singularity.
    def ode_system(x, Y):
        y, yp, ypp = Y  # y, y', y''
        
        # Denominator for the y'' term, ensuring it's not zero at x=0
        den = 3*x**2 + 2
        
        # The ODE: y''' = y^4 + (y')^4 - y''/(3x^2+2) + 0.5
        yppp = y**4 + yp**4 - ypp/den + 0.5
        return [yp, ypp, yppp]

    # Initial conditions
    y0 = [0.00, 3.00, 2.00]

    # Time span for integration. The singularity is expected to occur before x=1.
    t_span = [0, 1.0]

    # Numerically solve the ODE system. The solver will stop when it gets close
    # to the singularity, providing us with an estimate for x0 and y(x0).
    # Using a small max_step helps in accurately locating the singularity.
    sol = solve_ivp(ode_system, t_span, y0, method='RK45', max_step=1e-4, dense_output=True)

    # The singularity time x0 is the last time point reached by the solver.
    x0 = sol.t[-1]
    
    # The value A = y(x0) is the value of y at the singularity.
    A = sol.y[0, -1]

    # The coefficient B is derived analytically from the asymptotic analysis of y''' ~ (y')^4.
    B = (30)**(1/3)

    # Round the constants to two decimal places as requested.
    x0_r = round(x0, 2)
    A_r = round(A, 2)
    B_r = round(B, 2)

    print("The approximate analytical expression for y(x) in the large x regime (as x approaches the singularity) is:")
    final_expression = f"y(x) = {A_r:.2f} - {B_r:.2f}*({x0_r:.2f} - x)^(1/3)"
    print(final_expression)
    
    print("\nThe numbers in the final equation are:")
    print(f"A (the value of y(x) at the singularity) = {A_r:.2f}")
    print(f"B (the coefficient from asymptotic analysis) = {B_r:.2f}")
    print(f"x0 (the time of the singularity) = {x0_r:.2f}")

# Execute the function to get the final answer
solve_and_find_expression()
<<<y(x) = 2.51 - 3.11*(0.43 - x)^(1/3)>>>