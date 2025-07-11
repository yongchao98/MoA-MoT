import numpy as np

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using a Runge-Kutta (RK4) numerical method.
    """
    # 1. Define physical constants and parameters from the problem
    g = 9.8         # Acceleration due to gravity (m/s^2)
    d_c = 0.04      # Diameter of cardboard cylinder (m)
    r_c = d_c / 2   # Radius of cardboard cylinder (m)
    tau = 0.0005    # Thickness of a single paper layer (m)
    N = 100         # Number of wraps
    M_p = 0.200     # Mass of the paper (kg)
    M_c = 0.020     # Mass of the cardboard cylinder (kg)
    pi = np.pi

    # 2. Calculate derived properties of the roll
    # Total length of the paper (sum of circumferences of all layers)
    L = 2 * pi * (N * r_c + tau * N * (N + 1) / 2)
    # Linear mass density of the paper (kg/m)
    lambda_ = M_p / L
    # Square of the cylinder's radius
    r_c_sq = r_c**2

    # 3. Define functions for the roll's properties, which change as it unrolls.
    # These properties depend on 'y', the length of paper that has unrolled.
    
    def get_mass(y):
        """Calculates the total mass of the roll (cardboard + remaining paper)"""
        return M_c + lambda_ * (L - y)

    def get_radius_squared(y):
        """Calculates the square of the instantaneous outer radius of the roll"""
        # This is derived from Area = (L_rem * tau) = pi * (r_out^2 - r_c^2)
        return r_c_sq + (tau / pi) * (L - y)

    def get_moment_of_inertia(y):
        """Calculates the moment of inertia of the roll"""
        # I_c for a solid cylinder
        I_c = 0.5 * M_c * r_c_sq
        # I_p for a thick hollow cylinder
        m_p_remaining = lambda_ * (L - y)
        r_out_sq = get_radius_squared(y)
        I_p = 0.5 * m_p_remaining * (r_c_sq + r_out_sq)
        return I_c + I_p

    # 4. Define the function for acceleration using the energy conservation method.
    # The acceleration a=dv/dt is a function of position y and velocity v.
    
    def get_acceleration(y, v):
        """Calculates the acceleration a = dv/dt of the roll."""
        if y >= L:
            return 0.0

        # Effective mass of the system in the energy equation E = 0.5*M_eff*v^2 - U(y)
        m_val = get_mass(y)
        I_val = get_moment_of_inertia(y)
        r_sq_val = get_radius_squared(y)
        M_eff = m_val + I_val / r_sq_val
        
        # Derivative of potential energy U(y) w.r.t y is dU/dy = g*m(y)
        dU_dy = g * m_val
        
        # Derivative of effective mass w.r.t y (using quotient rule for I/r^2)
        dm_dy = -lambda_
        dr_sq_dy = -tau / pi
        dI_dy = -lambda_ * r_c_sq - lambda_ * (tau / pi) * (L - y)
        d_I_div_r_sq_dy = (dI_dy * r_sq_val - I_val * dr_sq_dy) / (r_sq_val**2)
        dM_eff_dy = dm_dy + d_I_div_r_sq_dy
        
        # From dE/dt = 0, we get a = (dU/dy - 0.5 * dM_eff/dy * v^2) / M_eff
        accel = (dU_dy - 0.5 * dM_eff_dy * v**2) / M_eff
        return accel
        
    # As requested, calculate and display the initial acceleration equation
    m0 = get_mass(0)
    I0 = get_moment_of_inertia(0)
    r0_sq = get_radius_squared(0)
    a0 = g * m0 / (m0 + I0 / r0_sq)
    print("The initial acceleration, a(0), is given by the equation:")
    print(f"a(0) = (g * m(0)) / (m(0) + I(0) / r(0)^2)")
    print(f"a(0) = ({g:.1f} * {m0:.3f}) / ({m0:.3f} + {I0:.6f} / {r0_sq:.6f}) = {a0:.2f} m/s^2\n")

    # 5. Implement the RK4 numerical solver
    # State vector S = [y, v], where y is position and v is velocity
    def dSdt(t, S):
        y, v = S
        return np.array([v, get_acceleration(y, v)])

    # Initial conditions
    y, v = 0.0, 0.0
    S = np.array([y, v])
    t = 0.0
    dt = 0.001  # Time step in seconds

    # Loop until the paper is fully unrolled (y >= L)
    while S[0] < L:
        S_prev = S
        
        k1 = dt * dSdt(t, S)
        k2 = dt * dSdt(t + 0.5*dt, S + 0.5*k1)
        k3 = dt * dSdt(t + 0.5*dt, S + 0.5*k2)
        k4 = dt * dSdt(t + dt, S + k3)
        
        S = S + (k1 + 2*k2 + 2*k3 + k4) / 6
        t = t + dt
        
    # 6. Interpolate to find the exact time when y = L
    y_prev = S_prev[0]
    y_now = S[0]
    t_prev = t - dt
    t_final = t_prev + dt * (L - y_prev) / (y_now - y_prev)
    
    print(f"The total length of the toilet paper is {L:.2f} meters.")
    print(f"The time it takes for the toilet paper to completely unroll is {t_final:.2f} seconds.")
    
    return t_final

# Run the solver and get the final answer
final_time = solve_toilet_paper_problem()
# Enclose final numerical answer in <<<>>>
# print(f"\n<<<{final_time:.2f}>>>") # This would be printed if not in a function
# Hack to get the final answer out for the format
global_final_answer = f"<<<{final_time:.2f}>>>"
print(f"\n{global_final_answer}")