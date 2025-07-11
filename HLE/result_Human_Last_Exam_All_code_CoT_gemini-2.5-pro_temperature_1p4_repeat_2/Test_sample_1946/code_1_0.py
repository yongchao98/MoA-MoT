import math

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a falling toilet paper roll to unravel.
    The problem is solved by numerically integrating the equations of motion using the RK4 method.
    """
    # 1. Define the constants of the problem
    g = 9.81  # Gravitational acceleration in m/s^2
    cardboard_diameter = 0.04  # m
    r_c = cardboard_diameter / 2  # Radius of the cardboard core in m
    d = 0.0005  # Thickness of the paper in m
    N = 100  # Number of wraps
    M_p = 0.2  # Mass of the paper in kg
    M_c = 0.02  # Mass of the cardboard core in kg

    # 2. Calculate derived properties of the roll
    # Calculate the initial outer radius of the full roll
    R_outer = r_c + N * d
    # Calculate the total length of the paper using the area method:
    # Area_side_view = L * thickness = pi * (R_outer^2 - r_c^2)
    L = math.pi / d * (R_outer**2 - r_c**2)

    # 3. Define functions to calculate roll properties as a function of unrolled length 'y'
    def get_acceleration(y):
        """Calculates the roll's acceleration 'a' for a given unrolled length 'y'."""
        
        # As paper unrolls, the outer radius 'r' decreases. We calculate r^2.
        # Based on area conservation: (L-y)*d = pi*(r^2 - r_c^2)
        r_sq = r_c**2 + (L - y) * d / math.pi
        
        # The mass of the paper remaining on the roll
        mass_paper_remaining = M_p * (1.0 - y / L)
        
        # The total mass of the falling roll (core + remaining paper)
        m_roll = M_c + mass_paper_remaining
        
        # The moment of inertia of the falling roll
        # I_roll = I_core + I_paper_remaining
        I_core = M_c * r_c**2
        # Paper is a thick hollow cylinder with inner radius r_c and outer radius r
        I_paper = 0.5 * mass_paper_remaining * (r_sq + r_c**2)
        I_roll = I_core + I_paper
        
        # The acceleration is derived from F=ma and tau=I*alpha
        # a = g * M / (M + I/r^2)
        acceleration = (g * m_roll) / (m_roll + I_roll / r_sq)
        
        return acceleration

    # 4. Implement the RK4 numerical integration method
    # Initial conditions
    t = 0.0
    y = 0.0  # Initial unrolled length (position)
    v = 0.0  # Initial velocity
    S = [y, v] # State vector [position, velocity]
    
    dt = 0.001  # Time step in seconds for the simulation

    # Loop until the entire roll has unraveled (y >= L)
    while S[0] < L:
        y_prev = S[0]
        t_prev = t

        # Standard RK4 algorithm
        # Derivatives are [dy/dt, dv/dt] = [v, a(y)]
        k1_v = S[1]
        k1_a = get_acceleration(S[0])
        
        k2_v = S[1] + 0.5 * dt * k1_a
        k2_a = get_acceleration(S[0] + 0.5 * dt * k1_v)
        
        k3_v = S[1] + 0.5 * dt * k2_a
        k3_a = get_acceleration(S[0] + 0.5 * dt * k2_v)
        
        k4_v = S[1] + dt * k3_a
        k4_a = get_acceleration(S[0] + dt * k3_v)

        # Update position and velocity
        S[0] += (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        S[1] += (dt / 6.0) * (k1_a + 2*k2_a + 2*k3_a + k4_a)
        
        t += dt

    # 5. Interpolate for a precise final time when y is exactly L
    y_curr = S[0]
    if y_curr > y_prev:
        t_final = t_prev + dt * (L - y_prev) / (y_curr - y_prev)
    else:
        t_final = t
        
    # Print the parameters and the final result
    print("--- Toilet Paper Fall Simulation ---")
    print("Initial Parameters:")
    print(f"  Gravity (g): {g} m/s^2")
    print(f"  Cardboard Core Diameter: {cardboard_diameter * 100} cm")
    print(f"  Cardboard Core Mass (M_c): {M_c} kg")
    print(f"  Total Paper Mass (M_p): {M_p} kg")
    print(f"  Paper Thickness (d): {d * 1000} mm")
    print(f"  Number of Wraps (N): {N}")
    
    print("\nCalculated Total Paper Length (L):")
    print(f"  L = (π/d) * (R_outer² - r_c²) = (π/{d}) * ({R_outer:.4f}² - {r_c}²) = {L:.2f} m")

    print("\nResult:")
    print("The time to unroll completely is found by numerically solving d²y/dt² = a(y), where a(y) is the acceleration.")
    print(f"Final calculated time: {t_final:.2f} seconds")
    
    # Return the numerical answer for the final tag
    return round(t_final, 2)

# Run the simulation and get the answer
final_time = solve_toilet_paper_fall()
# The final answer tag will be printed separately
# <<<final_time>>>