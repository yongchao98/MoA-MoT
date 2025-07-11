import scipy.integrate as spi
from fractions import Fraction

# The problem is to find the flux of the vector field F through the two yellow sides of a pyramid.
# The flux is given by the surface integral of F · dS.

# Vector Field F = (3x^3*y^2*z, 3x^2*y^3, z)

# The pyramid has its apex at (0,0,4) and base vertices at (+-1, +-1, 0).
# Let's define the first yellow face (S1) by vertices (0,0,4), (1,1,0), and (-1,1,0).
# This face lies on the plane z = 4 - 4y.
# The outward normal vector dS is (0, 4, 1) dx dy.
# The dot product F · dS = (3x^3*y^2*z, 3x^2*y^3, z) · (0, 4, 1)
#                        = 4 * (3x^2*y^3) + z
#                        = 12*x^2*y^3 + (4 - 4y)
def integrand_s1(x, y):
    """Integrand for the first yellow face."""
    return 12 * x**2 * y**3 + 4 - 4 * y

# The integration domain D1 for S1 is the projection onto the xy-plane,
# which is a triangle with vertices (0,0), (1,1), (-1,1).
# y ranges from 0 to 1, and for each y, x ranges from -y to y.
flux1, _ = spi.dblquad(integrand_s1, 0, 1, lambda y: -y, lambda y: y)

# The second yellow face (S3) is defined by vertices (0,0,4), (-1,-1,0), and (1,-1,0).
# This face lies on the plane z = 4 + 4y.
# The outward normal vector dS is (0, -4, 1) dx dy.
# The dot product F · dS = -4 * (3x^2*y^3) + z
#                        = -12*x^2*y^3 + (4 + 4y)
def integrand_s3(x, y):
    """Integrand for the second yellow face."""
    return -12 * x**2 * y**3 + 4 + 4 * y
    
# The integration domain D3 for S3 is the projection onto the xy-plane,
# which is a triangle with vertices (0,0), (-1,-1), (1,-1).
# y ranges from -1 to 0, and for each y, x ranges from y to -y.
flux3, _ = spi.dblquad(integrand_s3, -1, 0, lambda y: y, lambda y: -y)

# Due to symmetry, flux1 and flux3 are expected to be equal.
total_flux = flux1 + flux3

# For a more precise representation, we convert the results to fractions.
flux1_frac = Fraction(flux1).limit_denominator()
flux3_frac = Fraction(flux3).limit_denominator()
total_flux_frac = Fraction(total_flux).limit_denominator()

print("To find the energy flow, we calculate the flux through each of the two yellow sides and sum them.")
print(f"Flux through first yellow side = {flux1_frac.numerator}/{flux1_frac.denominator}")
print(f"Flux through second yellow side = {flux3_frac.numerator}/{flux3_frac.denominator}")
print("The final equation for the total energy flow is:")
print(f"{flux1_frac.numerator}/{flux1_frac.denominator} + {flux3_frac.numerator}/{flux3_frac.denominator} = {total_flux_frac.numerator}/{total_flux_frac.denominator}")
print(f"The numerical value of the total energy flow is approximately {total_flux:.4f}.")
