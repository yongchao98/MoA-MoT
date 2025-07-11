import sympy

def calculate_flux():
    """
    Calculates the flux of the vector field F through the two yellow sides of the pyramid.
    """
    # Define symbolic variables
    y, z = sympy.symbols('y z')

    # The problem specifies a square pyramid with its base on the z=0 plane centered at the origin,
    # with a base side of 2 and a height of 4.
    # The vertices of the base are (1,1,0), (1,-1,0), (-1,1,0), (-1,-1,0). The apex is (0,0,4).
    # The four side faces are on the planes:
    # S1 (right): 4x + z = 4
    # S2 (front): 4y + z = 4
    # S3 (left): -4x + z = 4
    # S4 (back): -4y + z = 4

    # The colors are "interspersed", meaning opposite faces have the same color.
    # This leads to an ambiguity: are the yellow faces (S1, S3) or (S2, S4)?
    # The vector field F = (3x^3*y^2*z, 3x^2*y^3, z) is not symmetric in x and y,
    # so the result will differ depending on the choice.
    # We proceed by choosing the yellow faces to be S1 and S3.

    # We will calculate the flux through S1 and multiply by 2, as the flux through S3 is identical due to symmetry.
    
    # For face S1 (4x + z = 4), we parametrize the surface using y and z.
    # x_s is x on the surface of S1
    x_s = 1 - z/4

    # The vector field F on the surface S1 is:
    # Fx = 3 * x_s**3 * y**2 * z
    # Fy = 3 * x_s**2 * y**3
    # Fz = z
    
    # The outward normal vector for a surface x=g(y,z) is (1, -dg/dy, -dg/dz).
    # Here, x = 1 - z/4. So g(y,z) = 1 - z/4.
    # The normal is (1, 0, 1/4).
    # dS is the surface element vector, which is (normal) * dy * dz.
    
    # F . dS = Fx * (1) + Fy * (0) + Fz * (1/4)
    # The integrand is the dot product of F and the normal vector component of dS.
    integrand = (3 * x_s**3 * y**2 * z) * 1 + (z) * (1/4)

    # The limits for y depend on z. The projection on the yz plane is a triangle.
    # For a given z, y ranges from -(1-z/4) to (1-z/4).
    y_limits = (y, -(1 - z/4), 1 - z/4)

    # Integrate with respect to y first
    inner_integral = sympy.integrate(integrand, y_limits)
    
    # Now integrate the result with respect to z from 0 to 4
    flux_s1 = sympy.integrate(inner_integral, (z, 0, 4))
    
    # By symmetry, the flux through the opposite face S3 is the same.
    flux_s3 = flux_s1
    
    # The total flux is the sum of the fluxes through the two yellow sides.
    total_flux = flux_s1 + flux_s3
    
    # Output the components of the final sum
    print(f"The flux through the first yellow side is: {flux_s1}")
    print(f"The flux through the second yellow side is: {flux_s3}")
    print(f"The total energy flow is the sum: {flux_s1} + {flux_s3} = {total_flux}")
    
calculate_flux()
<<<80/21>>>