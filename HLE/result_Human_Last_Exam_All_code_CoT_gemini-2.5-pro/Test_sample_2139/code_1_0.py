import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma

def solve_balloon_radius():
    """
    Solves the given differential equation for the radius of a spherical balloon y(t) at t=pi/4.
    """
    # 1. Define the ODE system dY/dt = f(t, Y), where Y = [y, y']
    def ode_system(t, Y):
        y, yp = Y[0], Y[1]

        s = np.sin(t)
        c = np.cos(t)
        
        # Avoid division by zero at t=0, though we start at t>0.
        if c == 0:
            # This case will not be hit in the integration range [epsilon, pi/4]
            return [0, 0]

        tan = s / c
        sec = 1.0 / c
        t4p1 = t**4 + 1
        
        # Coefficients from the differential equation: A*y'' + B*y' + C*y = F*y
        A = 4 * t4p1 * tan * sec
        B = 2 * (t**4 + 2 * tan * (t4p1 * tan + 8 * t**3) + 1) * sec
        C = 8 * t**2 * (t + 2 * tan * (t * tan + 3)) * sec
        F = t4p1 * np.sqrt(s)
        
        # Avoid division by zero for A if t=0.
        if A == 0:
            # This should only happen at t=0, which we are avoiding.
            return [yp, 0]

        # y'' = (F*y - B*y' - C*y) / A
        ypp = (F * y - B * yp - C * y) / A
        
        return [yp, ypp]

    # 2. Calculate the initial value y(0)
    gamma_2_3 = gamma(2.0/3.0)
    y0 = (128 * 3**(1.0/6.0) * gamma_2_3)**-1
    yp0 = 0.0
    
    # 3. Handle the singularity at t=0 by starting at a small epsilon
    # From small-t analysis: y(t) ≈ y0*(1 + t^(3/2)/6), y'(t) ≈ y0*t^(1/2)/4
    epsilon = 1e-9
    t_start = epsilon
    
    y_start = y0 * (1 + (1.0/6.0) * epsilon**1.5)
    yp_start = y0 * (1.0/4.0) * epsilon**0.5
    Y_start = [y_start, yp_start]

    # 4. Define the integration interval
    t_end = np.pi / 4.0
    t_span = [t_start, t_end]

    # 5. Solve the ODE
    # High precision is requested for a stable solution.
    sol = solve_ivp(
        ode_system, 
        t_span, 
        Y_start, 
        method='RK45', 
        dense_output=True, 
        rtol=1e-10, 
        atol=1e-12
    )

    # 6. Extract and print the final result
    # The solution at the last time point
    y_final = sol.y[0, -1]
    
    # The instruction "output each number in the final equation" is interpreted as
    # printing the key numbers in the final result statement.
    print(f"The initial radius is y(0) = {y0}")
    print(f"The final radius is y(pi/4) = {y_final}")
    
    # Final answer in the required format
    print(f"<<<{y_final}>>>")

solve_balloon_radius()