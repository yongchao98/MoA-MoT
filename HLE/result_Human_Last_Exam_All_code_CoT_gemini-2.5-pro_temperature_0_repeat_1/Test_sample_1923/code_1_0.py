import math

def solve_quadratic(a, b, c):
    """Solves ax^2 + bx + c = 0, returning the physically relevant (negative time) root."""
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        return None
    # We need the negative root for t_e, as t_e must be < 0 for an observation at T=0.
    # The term sqrt(b^2 - 4ac) is always greater than |b| in this physical setup.
    # So, to get a negative result, we need (-b - sqrt(...)) for positive 'a'.
    return (-b - math.sqrt(discriminant)) / (2 * a)

def calculate_shift():
    """
    Calculates the apparent shift of the center of gravity for a moving object
    based on the assumption in option C.
    """
    # --- Simulation Parameters ---
    # Speed of light
    c = 299792458.0  # m/s
    # Object's velocity (as a fraction of c)
    v_fraction = 0.5
    v = v_fraction * c  # m/s
    # Closest distance of approach (y-distance)
    d = 1.0e9  # meters
    # Length of the object (rod)
    L = 1.0e6  # meters
    epsilon = L / 2.0 # half-length

    print(f"--- Parameters ---")
    print(f"Observer at (0, 0)")
    print(f"Object moving at v = {v_fraction:.2f}c")
    print(f"Object's path is along the line y = {d:.1e} m")
    print(f"Object length = {L:.1e} m")
    print("-" * 20 + "\n")

    # --- Calculation for Leading Edge ---
    # Quadratic equation for emission time t_e: (c^2-v^2)t^2 - 2v*eps*t - (eps^2+d^2) = 0
    a = c**2 - v**2
    b_lead = -2 * v * epsilon
    c_const = -(epsilon**2 + d**2)
    t_e_lead = solve_quadratic(a, b_lead, c_const)
    t_prop_lead = -t_e_lead
    # Retarded position of the leading edge
    x_ret_lead = v * t_e_lead + epsilon
    y_ret = d # y-position is constant

    # --- Calculation for Trailing Edge ---
    # Quadratic equation for emission time t_e: (c^2-v^2)t^2 + 2v*eps*t - (eps^2+d^2) = 0
    b_trail = 2 * v * epsilon
    t_e_trail = solve_quadratic(a, b_trail, c_const)
    t_prop_trail = -t_e_trail
    # Retarded position of the trailing edge
    x_ret_trail = v * t_e_trail - epsilon

    # --- Apply Assumption from Option C ---
    # Force is proportional to 1 / t_prop
    F_lead = 1.0 / t_prop_lead
    F_trail = 1.0 / t_prop_trail

    # --- Calculate Center of Force ---
    # The center of force is the weighted average of the retarded positions,
    # weighted by the force magnitudes.
    center_of_force_x = (x_ret_lead * F_lead + x_ret_trail * F_trail) / (F_lead + F_trail)

    print("--- Results ---")
    print(f"Propagation time from leading edge:  {t_prop_lead:.8f} s")
    print(f"Propagation time from trailing edge: {t_prop_trail:.8f} s")
    print("(Note: Leading edge time is shorter)\n")

    print(f"Force from leading edge (prop. to 1/t):  {F_lead:.8e}")
    print(f"Force from trailing edge (prop. to 1/t): {F_trail:.8e}")
    print("(Note: Leading edge force is stronger)\n")

    print(f"Retarded X-position of leading edge:  {x_ret_lead:,.2f} m")
    print(f"Retarded X-position of trailing edge: {x_ret_trail:,.2f} m\n")

    print("--- Final Equation for Center of Force (X-coordinate) ---")
    print(f"x_cof = (x_ret_lead * F_lead + x_ret_trail * F_trail) / (F_lead + F_trail)")
    print(f"x_cof = ({x_ret_lead:.2f} * {F_lead:.3e} + {x_ret_trail:.2f} * {F_trail:.3e}) / ({F_lead:.3e} + {F_trail:.3e})")
    print(f"x_cof = {center_of_force_x:,.2f} m\n")

    print("--- Conclusion ---")
    if center_of_force_x > 0:
        print("The center of gravity appears shifted in the direction of motion (positive x).")
    else:
        print("The center of gravity does NOT appear shifted in the direction of motion.")

calculate_shift()