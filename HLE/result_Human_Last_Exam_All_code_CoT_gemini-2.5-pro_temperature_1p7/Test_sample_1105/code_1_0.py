import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and forces.

    This problem is solved using the principle of conservation of energy.
    An effective potential energy U_eff(r) = -F*r - G*m*M/r is defined.
    The total energy E = (1/2)M*v^2 + U_eff(r) is conserved.
    We find E from the initial conditions (l_0, v_0) and then solve for the
    maximum distance l_max where the velocity is zero, which leads to a
    quadratic equation: F*(l_max)^2 + E*l_max + G*m*M = 0.
    """

    # --- Constants and Initial Conditions ---
    # These values are chosen as a physically plausible example where a finite
    # maximum distance exists.

    # Gravitational constant in N * m^2 / kg^2
    G = 6.67430e-11
    # Mass of the asteroid A in kg
    m = 5.972e22  # (e.g., a large moon)
    # Mass of the spaceship B in kg
    M = 5.0e4     # (e.g., a large spacecraft)
    # Initial distance between A and B in meters
    l_0 = 8.0e6   # 8000 km
    # Initial speed of the spaceship B in m/s
    v_0 = 100
    # Additional constant force applied to B in Newtons
    F = 4000

    # --- Calculation Steps ---

    # 1. Calculate the total conserved energy (E_total) from initial conditions.
    try:
        term_kinetic = 0.5 * M * v_0**2
        term_potential_F = F * l_0
        term_potential_G = (G * m * M) / l_0
        E_total = term_kinetic - term_potential_F - term_potential_G
    except ZeroDivisionError:
        print("Error: Initial distance l_0 cannot be zero.")
        return

    # 2. Solve the quadratic equation for l_max:
    # a*x^2 + b*x + c = 0, where x = l_max
    a = F
    b = E_total
    c = G * m * M

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("The spaceship escapes to infinity (no real solution for a maximum distance).")
        print(f"This occurs because the total energy E ({E_total:.2e} J) is greater than the peak potential energy.")
        return

    # Calculate the two roots for l_max
    # We choose the larger root, which represents the outer turning point.
    l_max = (-b + math.sqrt(discriminant)) / (2 * a)

    # --- Final Output ---
    print("This script calculates the maximum distance l_max using the formula derived from conservation of energy.")
    print("The final formula to solve is: F*(l_max)^2 + E_total*l_max + G*m*M = 0\n")

    print("--- Input Values ---")
    print(f"G = {G} N m^2/kg^2")
    print(f"m = {m:.3e} kg")
    print(f"M = {M:.3e} kg")
    print(f"l_0 = {l_0:.3e} m")
    print(f"v_0 = {v_0} m/s")
    print(f"F = {F} N\n")

    print("--- Calculation ---")
    print("Total Energy E_total = (1/2)*M*v_0^2 - F*l_0 - G*m*M/l_0")
    print(f"E_total = (1/2)*{M:.3e}*{v_0}^2 - {F}*{l_0:.3e} - ({G}*{m:.3e}*{M:.3e})/{l_0:.3e}")
    print(f"E_total = {term_kinetic:.3e} - {term_potential_F:.3e} - {term_potential_G:.3e} = {E_total:.3e} Joules\n")
    
    print("Using the quadratic formula for l_max = [-E_total + sqrt(E_total^2 - 4*F*G*m*M)] / (2*F):")
    print(f"l_max = [-({E_total:.3e}) + sqrt(({E_total:.3e})^2 - 4*{F}*{c}) ] / (2*{F})")
    print(f"l_max = [{-E_total:.3e} + sqrt({discriminant:.3e}) ] / {2*F}")
    print(f"l_max = [{-E_total:.3e} + {math.sqrt(discriminant):.3e} ] / {2*F}\n")
    
    print("--- Result ---")
    print(f"The calculated maximum distance is: {l_max:.4f} meters")

if __name__ == '__main__':
    calculate_max_distance()