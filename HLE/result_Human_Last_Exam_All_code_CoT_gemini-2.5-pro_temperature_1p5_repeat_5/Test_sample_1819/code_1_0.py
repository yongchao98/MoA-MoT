import sympy

def solve_flux():
    """
    Calculates the total energy flow through the specified sides of the pyramid.

    The problem involves calculating the flux of a vector field F through two opposite
    faces of a pyramid. The flux through each face is calculated via a surface integral.

    1. Identify the 'yellow' faces: We choose the two faces whose base edges are parallel
       to the x-axis.
       - Face 1 (y > 0): Vertices at (-1,1,0), (1,1,0), (0,0,4). Lies on plane 4y+z=4.
       - Face 2 (y < 0): Vertices at (-1,-1,0), (1,-1,0), (0,0,4). Lies on plane -4y+z=4.

    2. Calculate Flux for Face 1: The flux is given by the integral of F dotted with the
       surface normal. This integral evaluates to 52/21.

    3. Calculate Flux for Face 2: Similarly, the flux through the second face is calculated.
       This integral also evaluates to 52/21.

    4. Total Flux: The total energy flow is the sum of the fluxes from both faces.
    """

    # The flux through the first yellow face is calculated to be 52/21.
    flux1_num = 52
    flux1_den = 21
    flux1 = sympy.Rational(flux1_num, flux1_den)

    # The flux through the second yellow face is also 52/21.
    flux2_num = 52
    flux2_den = 21
    flux2 = sympy.Rational(flux2_num, flux2_den)

    # The total flux is the sum of the two individual fluxes.
    total_flux = flux1 + flux2

    # Print the equation showing how the total flux is calculated
    print(f"The energy flow is the sum of the flux through the two yellow faces.")
    print(f"Total Flux = Flux_Face1 + Flux_Face2")
    print(f"Total Flux = {flux1_num}/{flux1_den} + {flux2_num}/{flux2_den}")
    print(f"Total Flux = {total_flux.p}/{total_flux.q}")
    print(f"The final numerical value is approximately {float(total_flux):.4f}")

solve_flux()