import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to fall and unroll completely.
    
    The problem is modeled using Newton's second law for a system with variable mass and
    moment of inertia. The resulting system of ordinary differential equations is solved
    numerically using the 4th-order Runge-Kutta (RK4) method.
    """
    # --- 1. Define Physical Parameters based on the problem statement ---
    
    # Diameter of the cardboard inner cylinder is 4 cm
    d_cyl_cm = 4.0 
    r_c = (d_cyl_cm / 100) / 2 # m, radius of the cylinder

    # Toilet paper thickness is 0.5 mm
    t_p_mm = 0.5
    t_p = t_p_mm / 1000 # m, thickness of paper

    # Number of wraps
    N = 100

    # Total mass of paper is 200 grams
    m_p_total_g = 200.0
    m_p_total = m_p_total_g / 1000 # kg
    
    # Mass of the cardboard cylinder is 20 grams
    m_c_g = 20.0
    m_c = m_c_g / 1000 # kg

    # Acceleration due to gravity
    g = 9.8

    # --- 2. Calculate Derived Constants ---

    # Total length of the paper, modeled as an Archimedean spiral.
    # L = integral(r d(theta)) = 2*pi*N*r_c + t_p*pi*N^2
    L_total = 2 * np.pi * N * r_c + t_p * np.pi * N**2

    # Mass per unit length of the paper
    lambda_p = m_p_total / L_total

    # --- 3. Define the System of Ordinary Differential Equations (ODEs) ---
    
    # The state of the system is a vector [y, v], where:
    # y = length of paper unrolled
    # v = downward velocity of the roll's center of mass
    def derivatives(state, t):
        y, v = state

        # To prevent math errors (e.g., division by zero) when y is very close
        # to L_total, we can stop the calculation if the roll is fully unrolled.
        if y >= L_total:
            return np.array([0.0, 0.0])

        # Calculate instantaneous properties of the falling roll
        
        # Mass of the falling part (core + remaining paper)
        M = m_c + lambda_p * (L_total - y)
        
        # Square of the outer radius of the remaining paper on the roll
        # This is derived from Area = pi*(R^2 - r_c^2) = Length_rem * thickness
        R_sq = r_c**2 + (t_p / np.pi) * (L_total - y)
        if R_sq <= 0: R_sq = 1e-12 # safety clamp to avoid sqrt of negative/zero
        
        # Moment of inertia of the falling part
        # I = I_core + I_paper
        # I_core for a thin hollow cylinder is m_c * r_c^2
        I_core = m_c * r_c**2
        # I_paper for a hollow cylinder is 1/2 * m_paper * (R^2 + r_c^2)
        m_paper_rem = lambda_p * (L_total - y)
        I_paper = 0.5 * m_paper_rem * (R_sq + r_c**2)
        I = I_core + I_paper
        
        # Calculate acceleration a = dv/dt using the full dynamics equation:
        # a = ( M*g + v^2*[dM/dy term] ) / ( M + I/R^2 )
        # The term v^2 * (lambda_p - I*t_p/(2*pi*R^4)) accounts for variable mass.
        numerator = M * g + v**2 * (lambda_p - (I * t_p) / (2 * np.pi * R_sq**2))
        denominator = M + I / R_sq
        
        a = numerator / denominator if denominator != 0 else 0.0

        # Return the derivatives [dy/dt, dv/dt]
        return np.array([v, a])

    # --- 4. Solve the ODEs using the RK4 method ---
    
    # Initial conditions
    y, v = 0.0, 0.0 
    t = 0.0
    dt = 0.001  # Time step in seconds

    while True:
        state = np.array([y, v])
        
        # RK4 calculation steps
        k1 = dt * derivatives(state, t)
        k2 = dt * derivatives(state + 0.5 * k1, t + 0.5 * dt)
        k3 = dt * derivatives(state + 0.5 * k2, t + 0.5 * dt)
        k4 = dt * derivatives(state + k3, t + dt)
        
        new_state = state + (k1 + 2*k2 + 2*k3 + k4) / 6.0
        
        # Check if the roll has finished unrolling in this time step
        if new_state[0] >= L_total:
            # Use linear interpolation to find the exact time it hit L_total
            y_prev = y
            y_next = new_state[0]
            fraction = (L_total - y_prev) / (y_next - y_prev)
            final_time = t + fraction * dt
            break
        
        # Update state and time for the next iteration
        y, v = new_state
        t += dt

    # --- 5. Print the results as requested ---
    
    print("This script solves for the time (t) it takes a toilet paper roll to unravel while falling.")
    print("It solves the system of differential equations dy/dt = v and dv/dt = a(y, v).")
    print("\n--- Input Parameters ---")
    print(f"Cylinder Diameter: {d_cyl_cm} cm")
    print(f"Paper Thickness: {t_p_mm} mm")
    print(f"Number of Wraps: {N}")
    print(f"Total Paper Mass: {m_p_total_g} g")
    print(f"Cardboard Core Mass: {m_c_g} g")
    print(f"Calculated Total Paper Length: {L_total:.2f} m")
    
    print("\n--- Final Result ---")
    # Final equation solved: Time = integral from 0 to L_total of (1/v) dy
    print(f"The equation for time is t = \u222B dy/v(y), from y=0 to y={L_total:.2f}.")
    print(f"The time for the toilet paper to completely unroll is {final_time:.2f} seconds.")
    
    return final_time

# Run the simulation and get the final answer
final_time_result = solve_toilet_paper_fall()
# The final answer is enclosed in <<<>>> as requested.
# print(f"\n<<<{final_time_result:.2f}>>>")