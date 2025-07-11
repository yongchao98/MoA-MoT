import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.

    The problem involves a ring of mass 'm' sliding on a horizontal rod,
    connected by a string of length 'L' to an object of mass 'M'. The
    tension is calculated at a specific angle 'theta' as the object falls.

    The final formula for tension T, derived from conservation of energy and
    Newton's laws, is:
    T = (M * m / (m + M * cos(theta)**2)) * (g * sin(theta) + L_omega_sq)
    where L_omega_sq is the term L*omega^2 found from energy conservation:
    L_omega_sq = 2 * g * sin(theta) / ( (m * sin(theta)**2) / (m + M) + cos(theta)**2 )
    """
    # --- Given parameters ---
    m = 1.0  # kg, mass of the ring
    M = 1.0  # kg, mass of the object
    g = 9.8  # m/s^2, acceleration due to gravity
    theta_deg = 60.0 # degrees

    # --- Step 1: Preliminaries ---
    # Convert angle to radians for use in trigonometric functions
    theta_rad = math.radians(theta_deg)
    s = math.sin(theta_rad)
    c = math.cos(theta_rad)

    # --- Step 2: Calculate the L*omega^2 term from energy conservation ---
    numerator_L_omega_sq = 2 * g * s
    denominator_L_omega_sq = (m * s**2) / (m + M) + c**2
    L_omega_sq = numerator_L_omega_sq / denominator_L_omega_sq

    # --- Step 3: Calculate the terms for the final tension equation ---
    # T = ( term_A / term_B ) * ( term_C + term_D )
    term_A = M * m
    term_B = m + M * c**2
    term_C = g * s
    term_D = L_omega_sq # This is L*omega^2

    # --- Step 4: Calculate the final tension ---
    T = (term_A / term_B) * (term_C + term_D)

    # --- Step 5: Output the results as requested ---
    print("This script calculates the tension in the string.")
    print(f"Given parameters: m={m} kg, M={M} kg, theta={theta_deg} degrees, g={g} m/s^2.")
    print("\nThe final tension is calculated using the formula: T = (A / B) * (C + D)")
    print("Where the components are:")
    print(f"  A = M * m = {term_A:.4f}")
    print(f"  B = m + M*cos^2(theta) = {term_B:.4f}")
    print(f"  C = g*sin(theta) = {term_C:.4f}")
    print(f"  D = L*omega^2 = {term_D:.4f}")
    print("\nPlugging in the numbers:")
    print(f"T = ({term_A:.4f} / {term_B:.4f}) * ({term_C:.4f} + {term_D:.4f})")
    print(f"The final tension in the string is: {T:.2f} Newtons")

# Execute the calculation
calculate_tension()
<<<28.52>>>