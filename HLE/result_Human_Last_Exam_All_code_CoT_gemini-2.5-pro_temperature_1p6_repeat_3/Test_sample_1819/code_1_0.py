from sympy import symbols, integrate, Rational

def calculate_flux():
    """
    This script calculates the energy flow through the yellow sides of a pyramid.

    The problem is to find the flux of the vector field F = (3x^3y^2z, 3x^2y^3, z)
    through two opposite sides of a square pyramid. The pyramid has its base on the
    z=0 plane, centered at the origin, with base side length 2, and height 4.

    The coloring is "interspersed", meaning opposite faces have the same color.
    The vector field F is not symmetric with respect to x and y, leading to two
    possible values for the flux, depending on which pair of faces is yellow:
    1. The two faces with bases on y=1 and y=-1.
    2. The two faces with bases on x=1 and x=-1.

    This script calculates the flux through the 'y-faces'. This choice is made
    as there is an ambiguity in the problem statement. The calculation relies on
    parameterizing the surfaces and performing a surface integral, which can be
    simplified using a change of variables (u,v -> s,t).
    """

    # The flux through a single 'y-face' (e.g., the one with its base on y=1)
    # can be broken down into two components after parameterization and integration.
    # Flux = Integral(Contribution from F_y) + Integral(Contribution from F_z)

    # Component from the F_y part of the vector field: 24 * integral( (u-v)^2 * (u+v)^3 )
    # after change of variables s=u+v, t=u-v, this integral becomes 8/7
    flux_comp_Fy = Rational(8, 7)

    # Component from the F_z part of the vector field: integral( 8 * (1-u-v) )
    # after change of variables, this integral becomes 4/3
    flux_comp_Fz = Rational(4, 3)

    # The flux through a single y-face is the sum of these two components.
    flux_single_y_face = flux_comp_Fy + flux_comp_Fz

    # Due to the symmetry of the setup and the field, the flux through the
    # opposite y-face (at y=-1) is identical.
    # The total energy flow through the two yellow y-faces is twice this value.
    total_flux_yellow_sides = 2 * flux_single_y_face

    print("The total energy flow is the sum of the fluxes through two opposite faces.")
    print("Let's assume the yellow faces are the 'y-faces'.")
    print(f"The flux through a single y-face is the sum of two components:")
    print(f"Component 1 (from F_y): {flux_comp_Fy.p}/{flux_comp_Fy.q}")
    print(f"Component 2 (from F_z): {flux_comp_Fz.p}/{flux_comp_Fz.q}")
    print(f"Flux per face = {flux_comp_Fy.p}/{flux_comp_Fy.q} + {flux_comp_Fz.p}/{flux_comp_Fz.q} = {flux_single_y_face.p}/{flux_single_y_face.q}")
    print(f"Total flux for two yellow faces = 2 * ({flux_single_y_face.p}/{flux_single_y_face.q}) = {total_flux_yellow_sides.p}/{total_flux_yellow_sides.q}")
    print(f"The final numerical value is {float(total_flux_yellow_sides)}")

calculate_flux()