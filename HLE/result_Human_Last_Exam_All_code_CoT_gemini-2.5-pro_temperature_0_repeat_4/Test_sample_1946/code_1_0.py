import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to unroll completely
    as it falls, using a 4th-order Runge-Kutta numerical method.
    """
    # 1. Define system parameters from the problem description
    g = 9.81  # Acceleration due to gravity (m/s^2)
    d_inner = 0.04  # Diameter of the inner cylinder (m)
    r_c = d_inner / 2.0 # Radius of the inner cylinder (m)
    t_paper = 0.0005  # Thickness of the paper (m)
    N = 100  # Number of wraps
    m_p = 0.2  # Mass of the paper (kg)
    m_c = 0.02  # Mass of the cardboard cylinder (kg)
    
    # 2. Calculate derived constants
    R_i = r_c + N * t_paper  # Initial outer radius
    # Total length of the paper, using a continuous model for consistency
    L = np.pi * (R_i**2 - r_c**2) / t_paper

    # 3. Define helper functions for y-dependent variables
    def get_R_sq(y):
        # Calculates the square of the outer radius as a function of unrolled length y
        if y >= L:
            return r_c**2
        return R_i**2 - (R_i**2 - r_c**2) * y / L

    def get_m_p_rem(y):
        # Calculates the remaining mass of the paper
        if y >= L:
            return 0.0
        return m_p * (1.0 - y / L)

    def get_M(y, m_p_rem):
        # Calculates the total mass of the falling roll
        return m_c + m_p_rem

    def get_I(y, R_sq, m_p_rem):
        # Calculates the moment of inertia of the falling roll
        I_c = m_c * r_c**2
        if m_p_rem <= 0:
            return I_c
        I_p = 0.5 * m_p_rem * (R_sq + r_c**2) # Thick hollow cylinder model
        return I_c + I_p

    def get_dX_dy(y, R_sq, m_p_rem):
        # Calculates the derivative of X = M + I/R^2 with respect to y
        if y >= L:
            return 0.0
        
        dR_sq_dy = -(R_i**2 - r_c**2) / L
        dm_p_rem_dy = -m_p / L
        
        # Using the chain rule and product rule for dX/dy
        # X = m_c + 1.5*m_p_rem + (m_c*r_c^2 + 0.5*m_p_rem*r_c^2) / R_sq
        term1 = 1.5 * dm_p_rem_dy
        d_inv_R_sq_dy = -1 / (R_sq**2) * dR_sq_dy
        term2 = (m_c * r_c**2) * d_inv_R_sq_dy
        term3_num = (dm_p_rem_dy * R_sq - m_p_rem * dR_sq_dy)
        term3 = 0.5 * r_c**2 * term3_num / (R_sq**2)
        
        return term1 + term2 + term3

    # 4. Define the function for the ODE system (dy/dt, dv/dt)
    def model(t, state):
        y, v = state
        
        if y >= L:
            return [0, 0]

        # Calculate all y-dependent quantities
        m_p_rem_y = get_m_p_rem(y)
        M_y = get_M(y, m_p_rem_y)
        R_sq_y = get_R_sq(y)
        I_y = get_I(y, R_sq_y, m_p_rem_y)
        
        # This is the effective inertial mass term from the yo-yo model
        X_y = M_y + I_y / R_sq_y
        
        # This term arises from the variable-mass nature of the system
        dX_dy_val = get_dX_dy(y, R_sq_y, m_p_rem_y)
        
        # Equation of motion for acceleration 'a'
        numerator = g * M_y - 0.5 * dX_dy_val * v**2
        denominator = X_y
        a = numerator / denominator
        
        return [v, a]

    # 5. Implement the 4th-order Runge-Kutta (RK4) solver
    y, v, t = 0.0, 0.0, 0.0
    dt = 0.001  # Time step in seconds

    state = np.array([y, v])

    while state[0] < L:
        k1 = np.array(model(t, state))
        k2 = np.array(model(t + 0.5*dt, state + 0.5*dt*k1))
        k3 = np.array(model(t + 0.5*dt, state + 0.5*dt*k2))
        k4 = np.array(model(t + dt, state + dt*k3))
        
        state_change = (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        
        # Store previous state for final interpolation
        y_prev = state[0]
        t_prev = t
        
        state += state_change
        t += dt

    # 6. Interpolate for a more accurate final time
    if (state[0] - y_prev) > 0:
        frac = (L - y_prev) / (state[0] - y_prev)
        final_t = t_prev + frac * dt
    else:
        final_t = t

    # 7. Print the results and final answer
    print("Problem Parameters:")
    print(f"Cylinder diameter: {d_inner*100:.0f} cm, Paper thickness: {t_paper*1000:.1f} mm, Wraps: {N}")
    print(f"Paper mass: {m_p*1000:.0f} g, Cylinder mass: {m_c*1000:.0f} g")
    print(f"Initial outer radius: {R_i*100:.1f} cm")
    print(f"Total paper length: {L:.2f} m")
    print("\nSimulation Results:")
    print(f"The final equation of motion is a complex ODE solved numerically.")
    print(f"The time it takes for the toilet paper to completely unroll is calculated to be:")
    print(f"Time = {final_t:.2f} seconds")

solve_toilet_paper_fall()
<<<3.79>>>