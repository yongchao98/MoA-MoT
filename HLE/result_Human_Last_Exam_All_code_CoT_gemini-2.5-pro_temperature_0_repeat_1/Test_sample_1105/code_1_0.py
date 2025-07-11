import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and applied forces.
    """
    # --- Input Parameters ---
    # Mass of the asteroid (kg)
    m = 6e24
    # Mass of the spaceship (kg)
    M = 5e4
    # Initial distance between asteroid and spaceship (m)
    l0 = 7e6
    # Initial speed of the spaceship (m/s)
    v0 = 1000
    # Additional force applied to the spaceship (N)
    F = 1e5
    # Universal gravitational constant (N*m^2/kg^2)
    G = 6.67430e-11

    print("--- Given Parameters ---")
    print(f"Mass of asteroid (m): {m:.2e} kg")
    print(f"Mass of spaceship (M): {M:.2e} kg")
    print(f"Initial distance (l0): {l0:.2e} m")
    print(f"Initial speed (v0): {v0:.2e} m/s")
    print(f"Applied force (F): {F:.2e} N")
    print(f"Gravitational constant (G): {G:.2e} N*m^2/kg^2\n")

    # The work-energy theorem gives us the equation:
    # F*(l_max - l0) - G*m*M*(1/l0 - 1/l_max) = -0.5*M*v0^2
    #
    # Rearranging this into a quadratic equation for l_max (x): a*x^2 + b*x + c = 0
    
    # Coefficient 'a'
    a = F
    
    # Coefficient 'b'
    # b = -(F*l0 + G*m*M/l0 - 0.5*M*v0^2)
    term1 = F * l0
    term2 = (G * m * M) / l0
    term3 = 0.5 * M * v0**2
    b = -(term1 + term2 - term3)

    # Coefficient 'c'
    c = G * m * M

    print("--- Solving the Quadratic Equation for l_max ---")
    print("Equation form: a*(l_max)^2 + b*(l_max) + c = 0")
    print(f"a = {a:.4e}")
    print(f"b = {b:.4e}")
    print(f"c = {c:.4e}\n")
    
    # The final equation with numerical values is:
    print("Final Equation:")
    print(f"({a:.4e}) * l_max^2 + ({b:.4e}) * l_max + ({c:.4e}) = 0\n")

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("The discriminant is negative. The spaceship has enough energy to escape to infinity.")
        print("Result: Maximum distance is infinite.")
        l_max = float('inf')
    else:
        # Solve for the two roots using the quadratic formula
        sqrt_discriminant = math.sqrt(discriminant)
        root1 = (-b - sqrt_discriminant) / (2 * a)
        root2 = (-b + sqrt_discriminant) / (2 * a)

        # The maximum distance is the larger of the two roots, as the spaceship moves away from l0.
        l_max = max(root1, root2)
        
        print("--- Result ---")
        print(f"The maximum distance (l_max) the spaceship can achieve is: {l_max:.4e} meters.")

    return l_max

if __name__ == "__main__":
    max_dist = calculate_max_distance()
    # The final answer is wrapped in <<<...>>>
    if max_dist != float('inf'):
        print(f"<<<{max_dist:.4e}>>>")
    else:
        print("<<<infinity>>>")
