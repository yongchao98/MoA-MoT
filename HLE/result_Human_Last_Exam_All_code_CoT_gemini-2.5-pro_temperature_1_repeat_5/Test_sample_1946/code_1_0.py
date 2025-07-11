import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a falling toilet paper roll to unravel.
    The problem is solved by numerically integrating the equations of motion
    using the fourth-order Runge-Kutta (RK4) method.
    """
    # 1. Define Constants and Parameters
    # -------------------------------------
    print("Physical Parameters Used in Calculation:")
    G = 9.81  # Acceleration due to gravity (m/s^2)
    R_C = 0.02  # Cardboard cylinder radius (m) (diameter is 4 cm)
    T_PAPER = 0.0005  # Toilet paper thickness (m) (0.5 mm)
    M_P = 0.2  # Mass of paper (kg) (200 g)
    M_C = 0.02  # Mass of cardboard cylinder (kg) (20 g)
    NUM_WRAPS = 100 # Number of times the paper is wrapped around the roll

    # Derived parameters
    # Initial outer radius of the full roll
    R_0 = R_C + NUM_WRAPS * T_PAPER
    # Total length of the paper, derived from the cross-sectional area
    # Area = pi*(R_0^2 - R_C^2) = L_TOTAL * T_PAPER
    L_TOTAL = np.pi * (R_0**2 - R_C**2) / T_PAPER

    print(f"  Gravity (g): {G} m/s^2")
    print(f"  Cardboard Radius (r_c): {R_C} m")
    print(f"  Paper Thickness (t): {T_PAPER} m")
    print(f"  Paper Mass (M_p): {M_P} kg")
    print(f"  Cardboard Mass (M_c): {M_C} kg")
    print(f"  Calculated Initial Outer Radius (R_0): {R_0:.4f} m")
    print(f"  Calculated Total Paper Length (L): {L_TOTAL:.4f} m\n")


    # 2. Formulate helper functions for the ODE
    # ------------------------------------------
    # These functions calculate the roll's properties as a function of
    # the unrolled length `y`.

    def get_radius_sq(y):
        """Calculates the square of the outer radius when length y is unrolled."""
        if y >= L_TOTAL: return R_C**2
        return R_C**2 + (L_TOTAL - y) * T_PAPER / np.pi

    def get_mass(y):
        """Calculates the mass of the remaining roll."""
        if y >= L_TOTAL: return M_C
        return M_C + M_P * (L_TOTAL - y) / L_TOTAL

    def get_moment_of_inertia(y):
        """Calculates the moment of inertia of the remaining roll."""
        if y >= L_TOTAL: return M_C * R_C**2
        
        # Moment of inertia of the cardboard (hollow cylinder, mass at R_C)
        I_c = M_C * R_C**2
        
        # Moment of inertia of remaining paper (hollow cylinder)
        # I_p = 0.5 * M_paper_rem * (r_outer^2 + r_inner^2)
        mass_paper_rem = M_P * (L_TOTAL - y) / L_TOTAL
        r_sq = get_radius_sq(y)
        I_p = 0.5 * mass_paper_rem * (r_sq + R_C**2)
        
        return I_c + I_p

    def get_acceleration(y, v):
        """
        Calculates the roll's acceleration for a given unrolled length y and velocity v.
        This function represents the differential equation dv/dt = a(y, v).
        """
        if y >= L_TOTAL: return 0.0

        m = get_mass(y)
        I = get_moment_of_inertia(y)
        r_sq = get_radius_sq(y)
        
        if r_sq <= 0: return 0.0

        # Full equation for acceleration including the v^2 term from dr/dt
        # a = (m*g - I*t*v^2 / (2*pi*r^4)) / (m + I/r^2)
        numerator = m * G - (I * T_PAPER * v**2) / (2 * np.pi * r_sq**2)
        denominator = m + I / r_sq
        
        if denominator == 0: return 0.0
        return numerator / denominator

    # 3. Numerical Simulation (RK4)
    # -----------------------------
    # Initial conditions
    t = 0.0
    y = 0.0  # Unrolled length
    v = 0.0  # Velocity

    # Simulation parameters
    dt = 0.001  # Time step in seconds

    while y < L_TOTAL:
        # Store previous state for final interpolation
        y_prev, t_prev = y, t

        # Standard RK4 algorithm for the system dy/dt = v, dv/dt = a(y,v)
        # Calculate slopes at intermediate points
        k1_y, k1_v = v, get_acceleration(y, v)
        
        k2_y = v + 0.5 * dt * k1_v
        k2_v = get_acceleration(y + 0.5 * dt * k1_y, v + 0.5 * dt * k1_v)
        
        k3_y = v + 0.5 * dt * k2_v
        k3_v = get_acceleration(y + 0.5 * dt * k2_y, v + 0.5 * dt * k2_v)

        k4_y = v + dt * k3_v
        k4_v = get_acceleration(y + dt * k3_y, v + dt * k3_v)

        # Update position, velocity, and time
        y += (dt / 6.0) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        v += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        t += dt

    # 4. Final Answer
    # ----------------
    # The loop stops when y > L_TOTAL. We linearly interpolate to find the
    # exact time when y was equal to L_TOTAL.
    if y > y_prev:
        final_time = t_prev + dt * (L_TOTAL - y_prev) / (y - y_prev)
    else: # Should not happen in this simulation
        final_time = t

    print(f"Final Result:")
    print(f"The time it takes for the toilet paper to completely unroll is {final_time:.2f} seconds.")
    
    # Return answer in the specified format
    return f"<<<{final_time:.2f}>>>"

# Run the simulation and print the final answer
final_answer = solve_toilet_paper_fall()
print(final_answer)
