import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and forces.
    """
    # Physical constants and parameters
    # Note: These values are chosen as an example for which a physical solution exists.
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    m = 1.0e12       # Mass of the asteroid A in kg
    M = 1.0e5        # Mass of the spaceship B in kg
    l_0 = 200.0      # Initial distance in meters
    v_0 = 0.1        # Initial speed in m/s
    F = 100.0        # Applied force in Newtons

    print("Given Parameters:")
    print(f"  Gravitational Constant (G): {G} m^3 kg^-1 s^-2")
    print(f"  Asteroid Mass (m): {m} kg")
    print(f"  Spaceship Mass (M): {M} kg")
    print(f"  Initial Distance (l_0): {l_0} m")
    print(f"  Initial Speed (v_0): {v_0} m/s")
    print(f"  Applied Force (F): {F} N")
    print("-" * 30)

    # The problem is solved by the work-energy theorem, leading to the equation:
    # F*l_max + GmM/l_max = F*l_0 + GmM/l_0 - 0.5*M*v_0^2
    # This can be written as a quadratic equation for l_max:
    # F*(l_max)^2 - C*(l_max) + GmM = 0
    # where C is the constant part of the first equation.

    # Calculate G*m*M for convenience
    GmM = G * m * M

    # Calculate the constant C
    try:
        C = F * l_0 + GmM / l_0 - 0.5 * M * v_0**2
    except ZeroDivisionError:
        print("Error: Initial distance l_0 cannot be zero.")
        return

    # Define coefficients of the quadratic equation a*x^2 + b*x + c = 0
    a = F
    b = -C
    c = GmM

    # Calculate the discriminant (D = b^2 - 4ac)
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("No real solution exists. The spaceship either escapes to infinity or the initial parameters do not allow for a turning point.")
        return

    # Calculate the two roots using the quadratic formula
    sqrt_discriminant = math.sqrt(discriminant)
    root1 = (-b - sqrt_discriminant) / (2 * a)
    root2 = (-b + sqrt_discriminant) / (2 * a)

    # The physical solution corresponds to the first turning point the spaceship encounters
    # as it moves away from l_0. This is the smallest root that is greater than l_0.
    # In a valid physical scenario for this problem, l_0 < root1 < root2.
    l_max = root1

    print("Derived Equation:")
    # We show the quadratic equation F*x^2 - C*x + GmM = 0 with numerical values
    # Let x = l_max
    print(f"  {a} * (l_max)^2 + ({b}) * l_max + ({c}) = 0")
    print("-" * 30)
    
    print("Final Calculation:")
    print(f"  The equation F*l_max + G*m*M/l_max = F*l_0 + G*m*M/l_0 - 1/2*M*v_0^2 becomes:")
    rhs = F * l_0 + GmM / l_0 - 0.5 * M * v_0**2
    print(f"  {F}*l_max + {GmM}/l_max = {F}*{l_0} + {GmM}/{l_0} - 1/2*{M}*{v_0**2}")
    print(f"  {F}*l_max + {GmM}/l_max = {rhs}")

    print("\nResult:")
    print(f"The maximum distance the spaceship can achieve is: {l_max:.4f} meters.")

calculate_max_distance()
<<<208.3060>>>