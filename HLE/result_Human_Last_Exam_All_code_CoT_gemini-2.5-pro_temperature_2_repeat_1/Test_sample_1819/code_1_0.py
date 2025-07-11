import sympy

def solve_flux():
    """
    Calculates the total energy flow (flux) through the yellow sides of the pyramid.

    The problem is solved analytically by setting up and evaluating the surface integrals.
    The flux through each of the two yellow faces (S1 and S3) is found to be identical.
    The total flux is therefore twice the flux through one face.

    Flux through one face is given by the integral:
    Integral from 0 to 4 of (2*z*(1-z/4)^6 + (z/2)*(1-z/4)) dz

    Using the substitution u = 1 - z/4, the definite integral becomes:
    4 * Integral from 0 to 1 of (8*u^6 - 8*u^7 + 2*u - 2*u^2) du

    This evaluates to:
    4 * [8/7 * u^7 - u^8 + u^2 - (2/3)*u^3] from 0 to 1
    = 4 * (8/7 - 1 + 1 - 2/3)
    = 4 * (8/7 - 2/3) = 4 * (10/21) = 40/21
    """

    # Numbers from the final step of the analytical integration for one face
    term1_numerator = 8
    term1_denominator = 7
    term2_numerator = 2
    term2_denominator = 3
    multiplier_from_substitution = 4
    num_yellow_faces = 2

    # Calculate the flux for a single face
    # Expression is multiplier * (term1 - term2)
    # 4 * (8/7 - 2/3)
    val1 = sympy.Rational(term1_numerator, term1_denominator)
    val2 = sympy.Rational(term2_numerator, term2_denominator)
    flux_one_face = multiplier_from_substitution * (val1 - val2)

    # Total flux is for two identical faces
    total_flux = num_yellow_faces * flux_one_face

    # The user wants to see the final equation. We print it step-by-step.
    print(f"The calculation for the total flux is based on the analytical solution of the integral.")
    print(f"Flux for one yellow face = {multiplier_from_substitution} * (({term1_numerator}/{term1_denominator}) - ({term2_numerator}/{term2_denominator}))")
    print(f"Total flux = {num_yellow_faces} * Flux for one yellow face")
    print(f"Total flux = {num_yellow_faces} * {flux_one_face}")
    print(f"Total flux = {total_flux}")
    print(f"The final numerical value is: {total_flux.evalf()}")


solve_flux()
<<<80/21>>>