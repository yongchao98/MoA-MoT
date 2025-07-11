from sympy import symbols, Matrix, integrate, Rational

def solve_flux():
    """
    Calculates the energy flow through the yellow sides of the pyramid.
    """
    # Define symbolic variables
    x, y, z = symbols('x y z')

    # Define the vector field F
    F = Matrix([3*x**3*y**2*z, 3*x**2*y**3, z])

    # We will calculate the flux through one of the yellow faces, S3, the "back" face.
    # Its vertices are (-1, 1, 0), (1, 1, 0), and (0, 0, 4).
    # The equation of the plane containing this face is y = 1 - z/4.
    
    # We parameterize the surface S3 using x and z as parameters.
    # r(x, z) = [x, 1 - z/4, z]
    # The parameter domain is the projection of the face onto the xz-plane,
    # which is a triangle with vertices (-1, 0), (1, 0), and (0, 4).
    # For integration, we let z run from 0 to 4, and for each z,
    # x runs from -(1 - z/4) to (1 - z/4).
    
    # Parameterize the surface
    r = Matrix([x, 1 - z/4, z])

    # Calculate the partial derivatives of r with respect to x and z
    dr_dx = r.diff(x)
    dr_dz = r.diff(z)

    # The normal vector dS is the cross product. We need the outward normal,
    # which for the back face must have a positive y-component.
    # dr_dx.cross(dr_dz) = [0, -1, -1/4], which is inward.
    # So we use dr_dz.cross(dr_dx).
    dS = dr_dz.cross(dr_dx) # Result: [0, 1, 1/4]

    # Substitute the parameterization into the vector field F
    F_on_S3 = F.subs(y, 1 - z/4)

    # Calculate the dot product F . dS for the integrand
    integrand = F_on_S3.dot(dS)

    # Define the integration limits for x
    w = 1 - z/4
    
    # Perform the integration to find the flux through one yellow face
    # First, integrate with respect to x
    flux_inner = integrate(integrand, (x, -w, w))
    
    # Then, integrate the result with respect to z
    flux_S3 = integrate(flux_inner, (z, 0, 4))

    # By symmetry, the flux through the other yellow face (S1, the "front" face)
    # is identical. The integrand F.n becomes the same when accounting for the
    # change in the sign of y and the direction of the normal vector.
    # Therefore, the total flux is twice the flux through S3.
    flux_S1 = flux_S3
    total_flux = flux_S1 + flux_S3

    # Print the final equation
    print(f"The flux through one yellow side is {flux_S3}.")
    print(f"The flux through the other yellow side is also {flux_S1}.")
    print(f"The total energy flow through the yellow sides is the sum of the fluxes.")
    print(f"Total Flow = {flux_S1} + {flux_S3} = {total_flux}")

solve_flux()