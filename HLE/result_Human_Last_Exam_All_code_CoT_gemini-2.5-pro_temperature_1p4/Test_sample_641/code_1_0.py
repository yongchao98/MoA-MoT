def solve_involutions():
    """
    Calculates the number of involutions in the projective special unitary group PSU(4, q).
    """
    # The problem specifies the finite field parameter q.
    q = 997

    # The number of involutions in PSU(4, q) for q odd is given by the formula:
    # N = (q^5 * (q^2 - q + 1) * (q^2 + 2*q - 1)) / 2
    # This is derived from counting preimages g in SU(4,q) where g^2 is central.

    # Calculate the components of the formula
    term1 = q**5
    term2 = q**2 - q + 1
    term3 = q**2 + 2 * q - 1

    # Calculate the total number of involutions.
    # We use integer division // as the result must be an integer.
    num_involutions = (term1 * term2 * term3) // 2

    # Output the result clearly, showing the formula and the calculated values.
    print(f"To find the number of involutions in PSU(4, {q}), we use the formula:")
    print(f"N = (q^5 * (q^2 - q + 1) * (q^2 + 2*q - 1)) / 2")
    print("\nSubstituting q = 997:")
    print(f"N = ({q}^5 * ({q}^2 - {q} + 1) * ({q}^2 + 2*{q} - 1)) / 2")
    print(f"N = ({term1} * {term2} * {term3}) / 2")
    print(f"\nThe number of involutions in PSU(4, 997) is:")
    print(num_involutions)

solve_involutions()