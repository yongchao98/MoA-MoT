import numpy as np
from scipy.integrate import solve_ivp

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to unroll completely while falling.

    The problem is modeled as an unrolling spool with variable mass, radius, and moment of inertia.
    A system of ordinary differential equations for the position and velocity of the roll is solved
    numerically using a Runge-Kutta method.
    """

    # 1. Define the physical constants from the problem description.
    # We will use g = 9.8 m/s^2 for the acceleration due to gravity.
    g = 9.8                   # Gravitational acceleration (m/s^2)
    d_inner_cm = 4.0            # Inner cylinder diameter (cm)
    t_paper_mm = 0.5            # Paper thickness (mm)
    m_paper_total_g = 200.0     # Total mass of paper (grams)
    m_cylinder_g = 20.0         # Mass of cardboard cylinder (grams)
    h_fall_m = 50.0             # Height of fall / total length of paper (m)

    # Convert units to SI (meters, kilograms, seconds)
    r_c = (d_inner_cm / 100) / 2.0  # Radius of cardboard cylinder (m)
    t_paper = t_paper_mm / 1000     # Paper thickness (m)
    m_paper_total = m_paper_total_g / 1000 # Total mass of paper (kg)
    m_c = m_cylinder_g / 1000       # Mass of cardboard cylinder (kg)
    L = h_fall_m                    # Total length of paper (m)

    # Derived constant: linear mass density of the paper
    rho_L = m_paper_total / L

    print("Solving for the time for the toilet paper to unroll, using a numerical simulation.")
    print("The model uses the following equation for acceleration a(y), where y is the unrolled length:")
    print("a(y) = g / (1 + I(y) / (M(y) * R(y)^2))\n")
    print("Parameters used in the final calculation:")
    print(f"  Gravitational acceleration (g): {g} m/s^2")
    print(f"  Inner cylinder radius (r_c): {r_c} m")
    print(f"  Paper thickness (t): {t_paper} m")
    print(f"  Total paper mass (m_paper_total): {m_paper_total} kg")
    print(f"  Cylinder mass (m_c): {m_c} kg")
    print(f"  Total paper length (L): {L} m\n")

    # 2. Define the system of differential equations: d(state)/dt = f(t, state)
    # The state is a vector [y, v], where y is position and v is velocity.
    def model(t, state):
        y, v = state

        # To avoid issues at the very end (y=L), we can handle it as a special case.
        if y >= L:
            return [0, 0] # Stop moving

        # Calculate M(y), R(y), and I(y)
        mass_paper_on_roll = rho_L * (L - y)
        M = m_c + mass_paper_on_roll

        # The cross-sectional area of the paper is conserved.
        # Area = (L-y)*t_paper = pi*(R^2 - r_c^2) => R^2 = r_c^2 + (L-y)*t_paper/pi
        R_sq = r_c**2 + (L - y) * t_paper / np.pi

        # Moment of inertia I = I_cylinder + I_paper
        # I_cylinder is for a thin ring: m_c * r_c^2
        # I_paper is for a thick hollow cylinder: 0.5 * m_paper_on_roll * (R^2 + r_c^2)
        I = m_c * r_c**2 + 0.5 * mass_paper_on_roll * (R_sq + r_c**2)

        # Acceleration a = g / (1 + I / (M * R^2))
        denominator = M * R_sq
        if denominator < 1e-12: # Avoid division by zero at the very end
             # At y=L, M=m_c, R_sq=r_c^2, I=m_c*r_c^2, so I/(M*R_sq)=1 -> a=g/2
             a = g / 2.0
        else:
             a = g / (1 + I / denominator)

        return [v, a]

    # 3. Set initial conditions and time span
    initial_conditions = [0.0, 0.0]  # Start at y=0 with v=0
    t_span = [0, 15] # Set a time span longer than expected

    # 4. Define an event to stop the integration when y = L
    def end_of_roll(t, state):
        return state[0] - L
    end_of_roll.terminal = True  # Stop the integration when this event occurs
    end_of_roll.direction = 1    # Trigger when event function goes from negative to positive

    # 5. Solve the ODE system
    solution = solve_ivp(
        fun=model,
        t_span=t_span,
        y0=initial_conditions,
        events=end_of_roll,
        method='RK45', # Standard Runge-Kutta method
        dense_output=True
    )

    # 6. Extract and print the result
    if solution.status == 1 and len(solution.t_events[0]) > 0:
        time_to_unroll = solution.t_events[0][0]
        print(f"The toilet paper takes {time_to_unroll:.2f} seconds to unroll completely.")
        # Final answer in the required format
        print(f"<<<{time_to_unroll:.2f}>>>")
    else:
        print("The simulation did not finish. The roll did not completely unroll in the given time.")
        print("<<<Error>>>")

solve_toilet_paper_fall()
