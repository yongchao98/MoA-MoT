import math
import numpy as np

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using a numerical ODE solver (RK4).
    """

    # Step 1: Define the physical constants from the problem description.
    # These are the "numbers" that define the equation of motion.
    g = 9.81  # Acceleration due to gravity (m/s^2)
    m_c = 0.020  # Mass of the cardboard cylinder (kg)
    m_p = 0.200  # Total mass of the paper (kg)
    d_c = 0.04   # Diameter of the cardboard cylinder (m)
    r_c = d_c / 2 # Radius of the cardboard cylinder (m)
    tau = 0.0005 # Thickness of the paper (m)
    N = 100      # Number of times the paper is wrapped around the roll

    # Step 2: Calculate derived properties of the roll.
    # Total initial mass of the roll
    M0 = m_c + m_p
    # Total length of the paper, calculated by summing the circumference of each wrap
    L_total = 2 * math.pi * (N * r_c + tau * N * (N + 1) / 2)
    # Initial outer radius of the full roll
    R0 = r_c + N * tau
    # Linear mass density of the paper
    lambda_L = m_p / L_total

    # The equation being solved is a second-order ODE: d^2l/dt^2 = a(l, dl/dt)
    # The parameters defining this specific problem instance are:
    print("--- System Parameters ---")
    print(f"Gravity (g): {g} m/s^2")
    print(f"Cardboard Mass (m_c): {m_c} kg")
    print(f"Paper Mass (m_p): {m_p} kg")
    print(f"Cardboard Radius (r_c): {r_c} m")
    print(f"Paper Thickness (tau): {tau} m")
    print(f"Number of Wraps (N): {N}")
    print(f"Calculated Total Length (L_total): {L_total:.4f} m")
    print(f"Calculated Initial Radius (R0): {R0:.4f} m")
    print("--------------------------\n")
    
    # Step 3: Define the function that calculates the derivatives [dl/dt, dv/dt]
    # This function represents the ODE: [v, a(l, v)]
    def get_derivatives(state, t):
        l, v = state # l: unrolled length, v: velocity

        if l >= L_total:
            return np.array([0.0, 0.0])

        # Calculate instantaneous properties of the roll
        L_rem = L_total - l
        m_p_rem = m_p * (L_rem / L_total)
        M = m_c + m_p_rem # Total current mass of the roll
        
        # Avoid division by zero or sqrt of negative, though unlikely with l < L_total
        R_sq_arg = L_rem / L_total
        if R_sq_arg < 0: R_sq_arg = 0
        R_sq = r_c**2 + (R0**2 - r_c**2) * R_sq_arg
        R = math.sqrt(R_sq)
        
        # Moment of inertia for the cardboard (thin hoop: I = mr^2)
        I_c = m_c * r_c**2
        # Moment of inertia for the paper (thick hollow cylinder: I = 1/2*m*(R_outer^2 + R_inner^2))
        I_p = 0.5 * m_p_rem * (R_sq + r_c**2)
        I = I_c + I_p # Total current moment of inertia

        # Calculate time derivatives of M, R, and I needed for the acceleration formula
        M_dot = -lambda_L * v
        
        if R == 0:
            R_dot = 0.0
        else:
            R_dot = -((R0**2 - r_c**2) * v) / (2 * R * L_total)

        # dI/dt = d/dt [I_c + 0.5 * m_p_rem * (R^2 + r_c^2)]
        m_p_rem_dot = -lambda_L * v
        R_sq_dot = 2 * R * R_dot
        I_dot = 0.5 * (m_p_rem_dot * (R_sq + r_c**2) + m_p_rem * R_sq_dot)
        
        # This is the full equation for acceleration, derived from Newton's laws
        # a = (g*M*R^2 - v*(dM/dt*R^2 + dI/dt - I/R*dR/dt)) / (M*R^2 + I)
        v_term = M_dot * R_sq + I_dot - I * R_dot / R
        numerator = g * M * R_sq - v * v_term
        denominator = M * R_sq + I

        if denominator == 0:
            a = 0.0
        else:
            a = numerator / denominator

        return np.array([v, a])

    # Step 4: Perform numerical integration using RK4
    t = 0.0
    # state is a numpy array: [unrolled_length, velocity]
    state = np.array([0.0, 0.0])
    dt = 0.001  # Time step for the simulation

    while state[0] < L_total:
        state_old = state
        t_old = t

        k1 = get_derivatives(state, t)
        k2 = get_derivatives(state + 0.5 * dt * k1, t + 0.5 * dt)
        k3 = get_derivatives(state + 0.5 * dt * k2, t + 0.5 * dt)
        k4 = get_derivatives(state + dt * k3, t + dt)

        state = state + (k1 + 2*k2 + 2*k3 + k4) / 6.0 * dt
        t = t + dt

    # Step 5: Interpolate to find the exact time the roll finishes
    l_prev = state_old[0]
    l_curr = state[0]
    
    if l_curr > l_prev:
        # t_final = t_old + dt * (fraction of last step needed)
        t_final = t_old + dt * (L_total - l_prev) / (l_curr - l_prev)
    else:
        t_final = t

    print("--- Final Result ---")
    print(f"The final equation is the time 't' when the unrolled length l(t) = L_total.")
    print(f"Solving l(t) = {L_total:.4f} m gives...")
    print(f"\nTime to unroll completely: {t_final:.2f} s")
    return t_final

final_time = solve_toilet_paper_problem()
<<<3.07>>>