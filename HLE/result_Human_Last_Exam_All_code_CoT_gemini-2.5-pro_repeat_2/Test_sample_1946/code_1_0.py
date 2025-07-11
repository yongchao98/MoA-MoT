import numpy as np
import math

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using a numerical simulation.
    """
    # --- 1. Define Constants and Parameters ---
    g = 9.8  # Gravitational acceleration (m/s^2)
    
    # Cardboard core parameters
    d_c = 4.0  # cm
    r_c = d_c / 2 / 100  # Core radius (m)
    m_c = 20.0 / 1000  # Core mass (kg)
    
    # Paper parameters
    t_paper = 0.5 / 1000  # Paper thickness (m)
    N_wraps = 100
    m_p_total = 200.0 / 1000  # Total paper mass (kg)
    
    # Initial conditions
    R_0 = r_c + N_wraps * t_paper  # Initial outer radius (m)
    v_0 = 0.0  # Initial velocity (m/s)
    y_0 = 0.0  # Initial position (m)

    # Derived constants for calculations
    # Moment of inertia of the core (thin cylindrical shell model: I = m*r^2)
    I_c = m_c * r_c**2
    # This term is constant in mass and inertia calculations
    R0_sq_minus_rc_sq = R_0**2 - r_c**2

    # --- 2. Print Model Setup ---
    print("--- Physics Model and Parameters ---")
    print(f"Gravitational acceleration (g): {g} m/s^2")
    print(f"Cardboard core radius (r_c): {r_c:.3f} m")
    print(f"Cardboard core mass (m_c): {m_c:.3f} kg")
    print(f"Paper thickness (t): {t_paper:.4f} m")
    print(f"Total paper mass (m_p_total): {m_p_total:.3f} kg")
    print(f"Initial roll radius (R_0): {R_0:.3f} m")
    
    print("\n--- Equations of Motion ---")
    print("The system is described by a set of differential equations solved numerically.")
    print("State Vector S = [position, velocity, radius]")
    print("dS/dt = [v, a(R), dR/dt]")
    print("\n1. Downward Acceleration `a(R)`:")
    print("   a(R) = M(R) * g / (M(R) + I(R) / R^2)")
    print("   Where:")
    print(f"     M(R) = m_c + m_p_total * (R^2 - r_c^2) / ({R0_sq_minus_rc_sq:.5f})  (Instantaneous Mass)")
    print(f"     I(R) = {I_c:.2e} + 0.5 * m_p_total * (R^4 - r_c^4) / ({R0_sq_minus_rc_sq:.5f})  (Instantaneous Moment of Inertia)")
    print("\n2. Rate of Radius Change `dR/dt`:")
    print(f"   dR/dt = -v * {t_paper:.4f} / (2 * pi * R)")
    
    # --- 3. Define Functions for the ODE System ---
    def get_acceleration(R):
        """Calculates instantaneous acceleration based on radius R."""
        if R <= r_c:
            return 0
        
        # Instantaneous mass of the paper
        m_p_R = m_p_total * (R**2 - r_c**2) / R0_sq_minus_rc_sq
        # Instantaneous total mass of the roll
        M_R = m_c + m_p_R
        
        # Instantaneous moment of inertia of the paper
        I_p_R = 0.5 * m_p_total * (R**4 - r_c**4) / R0_sq_minus_rc_sq
        # Instantaneous total moment of inertia
        I_R = I_c + I_p_R
        
        # Acceleration
        denominator = M_R + I_R / R**2
        if denominator == 0:
            return 0
        a = (M_R * g) / denominator
        return a

    def derivatives(state):
        """Calculates the derivatives [dy/dt, dv/dt, dR/dt] for the RK4 solver."""
        y, v, R = state
        if R <= r_c:
            return np.array([0.0, 0.0, 0.0])
            
        a = get_acceleration(R)
        dR_dt = -v * t_paper / (2 * math.pi * R)
        
        return np.array([v, a, dR_dt])

    # --- 4. Run the RK4 Simulation ---
    # Simulation parameters
    dt = 0.001  # Time step (s)
    
    # Initialize variables
    time = 0.0
    # State vector: [position, velocity, radius]
    state = np.array([y_0, v_0, R_0])
    
    state_prev = state
    time_prev = time

    while state[2] > r_c:
        # Store previous state for final interpolation
        state_prev = np.copy(state)
        time_prev = time
        
        # RK4 step
        k1 = derivatives(state)
        k2 = derivatives(state + 0.5 * dt * k1)
        k3 = derivatives(state + 0.5 * dt * k2)
        k4 = derivatives(state + dt * k3)
        
        state += (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        time += dt

    # --- 5. Interpolate for the exact final time ---
    R_prev = state_prev[2]
    R_curr = state[2]
    
    # t_final = t_prev + dt * (target - y_prev) / (y_curr - y_prev)
    if (R_curr - R_prev) != 0:
        t_final = time_prev + dt * (r_c - R_prev) / (R_curr - R_prev)
    else:
        t_final = time

    print("\n--- Simulation Result ---")
    print(f"The simulation stopped when the roll radius reached the core radius.")
    print(f"Time to unroll completely: {t_final:.2f} seconds")
    return t_final

if __name__ == '__main__':
    final_time = solve_toilet_paper_problem()
    print(f"\n<<< {final_time:.2f} >>>")
