import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    given an initial velocity and a constant applied force.
    """
    # --- Given Parameters (using example SI unit values) ---
    # Asteroid mass (kg)
    m = 5.972e22
    # Spaceship mass (kg)
    M = 1.0e5
    # Initial distance (m)
    l_0 = 1.0e6
    # Initial speed (m/s)
    v_0 = 100.0
    # Applied force (N)
    F = 1000.0
    # Gravitational constant (N m^2 / kg^2)
    G = 6.674e-11

    # --- Physics Formulation ---
    # The problem is solved using the work-energy theorem:
    # Change in Mechanical Energy = Work by Non-Conservative Forces
    # E_final - E_initial = W_F
    #
    # At l_max, the final velocity is 0.
    # (-G*m*M / l_max) - (0.5*M*v_0**2 - G*m*M / l_0) = F * (l_max - l_0)
    #
    # This can be rearranged into a quadratic equation for l_max:
    # F * l_max^2 - (F*l_0 + G*m*M/l_0 - 0.5*M*v_0**2) * l_max + G*m*M = 0
    # This is in the form: a*x^2 + b*x + c = 0, where x = l_max.

    print("--- Input Parameters ---")
    print(f"Asteroid Mass (m): {m:.3e} kg")
    print(f"Spaceship Mass (M): {M:.3e} kg")
    print(f"Initial Distance (l_0): {l_0:.3e} m")
    print(f"Initial Speed (v_0): {v_0:.1f} m/s")
    print(f"Applied Force (F): {F:.1f} N\n")

    # --- Calculation Steps ---
    # 1. Define the coefficients of the quadratic equation: a*x^2 + b*x + c = 0
    a = F
    
    # Calculate the term C for the 'b' coefficient
    # C = F*l_0 + G*m*M/l_0 - 0.5*M*v_0**2
    C = F * l_0 + (G * m * M) / l_0 - 0.5 * M * v_0**2
    b = -C
    
    c = G * m * M

    print("--- Solving the Quadratic Equation for l_max ---")
    print("Equation form: a*(l_max)^2 + b*(l_max) + c = 0")
    print(f"a = F = {a:.3e}")
    print(f"b = -(F*l_0 + G*m*M/l_0 - 0.5*M*v_0^2) = {b:.3e}")
    print(f"c = G*m*M = {c:.3e}\n")

    # 2. Calculate the discriminant (D = b^2 - 4*a*c)
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("Result: The discriminant is negative.")
        print("This means the spaceship has enough energy to escape to infinity.")
        print("There is no finite maximum distance.")
        return

    # 3. Solve for l_max using the quadratic formula: x = (-b ± sqrt(D)) / 2a
    # We take the '+' sign for the larger root, which corresponds to the maximum distance.
    l_max = (-b + math.sqrt(discriminant)) / (2 * a)

    # --- Final Output ---
    print("The maximum distance l_max is found using the quadratic formula:")
    print("l_max = (-b + sqrt(b^2 - 4*a*c)) / (2*a)\n")
    print("--- Plugging in the numbers ---")
    print(f"l_max = ( -({b:.3e}) + sqrt( ({b:.3e})^2 - 4 * ({a:.3e}) * ({c:.3e}) ) ) / ( 2 * {a:.1f} )")
    print(f"l_max = ( {(-b):.3e} + sqrt( {discriminant:.3e} ) ) / ( {2*a:.1f} )")
    print(f"l_max = ( {(-b):.3e} + {math.sqrt(discriminant):.3e} ) / ( {2*a:.1f} )")
    print(f"l_max = ( {(-b + math.sqrt(discriminant)):.3e} ) / ( {2*a:.1f} )\n")

    print(f"Calculated Maximum Distance (l_max): {l_max:.4e} meters")
    
    # Return the final value for the answer tag
    return l_max

if __name__ == '__main__':
    final_l_max = calculate_max_distance()
    if final_l_max:
        print(f"\n<<<__{final_l_max:.4e}__>>>")
