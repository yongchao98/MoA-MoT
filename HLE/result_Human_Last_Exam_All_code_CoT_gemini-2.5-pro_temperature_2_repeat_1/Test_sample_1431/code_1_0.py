def solve_rope_challenge():
    """
    Calculates the force F required to lift a rope of mass m and length l,
    such that its final velocity is v when fully lifted.
    """
    # --- Given Parameters ---
    # You can change these values to solve for a different scenario
    m = 10.0   # mass of the rope in kg
    l = 20.0   # length of the rope in m
    v = 5.0    # final speed of the rope in m/s
    g = 9.8    # acceleration due to gravity in m/s^2

    # --- Physics Explanation ---
    print("To solve the monk's challenge, we use the Work-Energy Theorem.")
    print("The theorem states: Work_done_by_Force_F + Work_done_by_Gravity = Change_in_Kinetic_Energy")
    print(f"This translates to the equation: (F * l) - (m * g * l / 2) = (1/2) * m * v^2\n")

    print("First, let's calculate the energy changes:")
    
    # Calculate Change in Potential Energy
    delta_pe = m * g * (l / 2)
    print(f"1. Change in Potential Energy (ΔPE = m*g*l/2): {m} * {g} * ({l}/2) = {delta_pe:.2f} J")
    
    # Calculate Change in Kinetic Energy
    delta_ke = 0.5 * m * v**2
    print(f"2. Change in Kinetic Energy (ΔKE = 1/2*m*v²): 0.5 * {m} * {v}² = {delta_ke:.2f} J\n")
    
    # The work done by the applied force F must equal the sum of these energy gains.
    print("The work done by the force F must equal the sum of ΔPE and ΔKE.")
    print(f"Work_F = {delta_pe:.2f} J + {delta_ke:.2f} J = {delta_pe + delta_ke:.2f} J")
    print(f"Since Work_F = F * l, we can solve for F: F = Work_F / l\n")

    # --- Final Calculation ---
    print("The final equation for the force F is: F = (m/2) * (v²/l + g)")
    print("Plugging in the numbers:")

    # Calculate F
    force = (m / 2) * (v**2 / l + g)

    # Print the equation with all the numbers
    print(f"F = ({m} / 2) * ({v}² / {l} + {g})")
    print(f"F = {m / 2} * ({(v**2 / l):.2f} + {g})")
    print(f"F = {m / 2} * {(v**2 / l + g):.2f}")
    print(f"F = {force:.2f} N\n")

    # --- The Answer ---
    print("The precise mystical force required is:")
    print(force)
    
solve_rope_challenge()