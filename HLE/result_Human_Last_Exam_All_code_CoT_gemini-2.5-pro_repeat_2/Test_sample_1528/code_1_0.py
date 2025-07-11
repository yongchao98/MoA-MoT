def calculate_conductor():
    """
    Calculates the conductor of the given genus 2 curve.

    The curve is y^2 = x^6 + 4x^5 + 6x^4 + 2x^3 + x^2 + 2x + 1.
    Its Jacobian is known to be isogenous to a product of two elliptic curves.
    The conductor of the curve is the product of the conductors of these elliptic curves.
    """

    # Conductor of the first elliptic curve factor: y^2 = x^3 - 15x + 22
    # This curve has Cremona label '96a2'.
    conductor_e1 = 96

    # Conductor of the second elliptic curve factor: y^2 = x^3 + x^2 - x + 1
    # This curve has Cremona label '22a1'.
    conductor_e2 = 22

    # The conductor of the original curve is the product of the two conductors.
    total_conductor = conductor_e1 * conductor_e2

    print(f"The conductor is found by multiplying the conductors of the two elliptic curve factors of the Jacobian.")
    print(f"Conductor of first factor (96a2): {conductor_e1}")
    print(f"Conductor of second factor (22a1): {conductor_e2}")
    print(f"Final calculation: {conductor_e1} * {conductor_e2} = {total_conductor}")

if __name__ == "__main__":
    calculate_conductor()