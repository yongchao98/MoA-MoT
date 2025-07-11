import math

def solve_spaceship_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    given an initial velocity and an additional constant force.
    """
    # Define the physical constants and variables.
    # These values are chosen to demonstrate a scenario where a finite
    # maximum distance exists.
    # Gravitational constant in m^3 kg^-1 s^-2
    G = 6.67430e-11
    # Mass of asteroid A (kg) - using a value comparable to a large moon
    m = 7.342e22
    # Mass of spaceship B (kg)
    M = 5e4
    # Initial distance (m) - starting from a low orbit
    l0 = 1.8e6
    # Initial speed (m/s)
    v0 = 1600.0
    # Additional constant force (N) - must be less than gravity at l0 for this problem
    F = 1.5e5

    print("Given parameters:")
    print(f"Gravitational constant G = {G:.5e} N m^2/kg^2")
    print(f"Mass of asteroid m = {m:.3e} kg")
    print(f"Mass of spaceship M = {M:.3e} kg")
    print(f"Initial distance l0 = {l0:.3e} m")
    print(f"Initial speed v0 = {v0:.1f} m/s")
    print(f"Additional force F = {F:.3e} N\n")

    GMm = G * M * m
    
    # Check for escape conditions. A finite maximum distance only exists if the
    # spaceship doesn't have enough thrust or velocity to escape.
    
    # Condition 1: Net force must be retarding. F < G*M*m/l^2
    l_crit = math.sqrt(GMm / F)
    if l0 > l_crit:
        print("Escape Condition Met:")
        print(f"The critical distance where thrust equals gravity is {l_crit:.3e} m.")
        print(f"Since l0 ({l0:.3e} m) > l_crit, the net force is always outwards.")
        print("The spaceship will escape to infinity. No finite maximum distance exists.")
        return

    # From the work-energy theorem, we derive the quadratic equation:
    # F * l_max^2 - K * l_max + G*M*m = 0
    # where K is a constant determined by the initial conditions.
    K = F * l0 + GMm / l0 - 0.5 * M * v0**2

    # The discriminant of this equation must be non-negative for real roots to exist.
    discriminant = K**2 - 4 * F * GMm

    if discriminant < 0:
        print("Escape Condition Met:")
        print("The initial velocity is too high. The spaceship has enough energy to")
        print("overcome the potential barrier and will escape to infinity.")
        print("No finite maximum distance exists.")
        return

    # Solve the quadratic equation for l_max. The maximum distance reached is
    # the smaller of the two roots, representing the first turning point.
    sqrt_discriminant = math.sqrt(discriminant)
    l_max = (K - sqrt_discriminant) / (2 * F)

    print("The final equation for the maximum distance l_max is the solution to:")
    print(f"{F:.4e} * l_max^2 - ({K:.4e}) * l_max + {GMm:.4e} = 0")
    print("\nUsing the quadratic formula, we find:")
    final_equation_str = (
        f"l_max = [ {K:.4e} - sqrt(({K:.4e})^2 - 4 * {F:.4e} * {GMm:.4e}) ] / (2 * {F:.4e})"
    )
    print(final_equation_str)

    print(f"\nl_max = [ {K:.4e} - sqrt({discriminant:.4e}) ] / {2 * F:.4e}")
    print(f"l_max = [ {K:.4e} - {sqrt_discriminant:.4e} ] / {2 * F:.4e}")
    print(f"l_max = {l_max:.4e} meters")
    print("\n-----------------------------------------------------")
    print(f"The maximum distance the spaceship can achieve is {l_max:.4e} meters.")
    print("-----------------------------------------------------")


# Execute the function to solve the problem
solve_spaceship_distance()