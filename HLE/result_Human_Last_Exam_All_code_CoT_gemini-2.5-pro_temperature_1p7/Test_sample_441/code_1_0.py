import numpy as np

def get_R_phi():
    """
    Calculates the parameters R and phi for the deferent-epicycle model of an object
    moving around a square.

    The derivation is based on Fourier analysis of the object's path.
    """

    # The deferent corresponds to the Fourier component with the largest magnitude,
    # and the epicycle to the one with the second-largest magnitude.
    # From the Fourier analysis, the non-zero coefficients are for n = 1, -3, 5, -7, ...
    # The magnitude of the coefficients is proportional to 1/n^2.
    # Thus, the largest coefficient is for n=1 (deferent) and the second largest for n=-3 (epicycle).

    n_def = 1
    n_epi = -3

    # The magnitude of the Fourier coefficient c_n is given by:
    # |c_n| = (8 * sqrt(2)) / (n^2 * pi^2)
    # We define a function for the scaled magnitude, omitting common constant factors.
    def c_mag_scaled(n):
        return 1 / (n**2)

    # The radius of the deferent (R_def) and epicycle (r_epi) are proportional to
    # the magnitude of their respective Fourier coefficients.
    R_def_scaled = c_mag_scaled(n_def)
    r_epi_scaled = c_mag_scaled(n_epi)

    # R is the ratio of the radius of the deferent to the radius of the epicycle.
    R = R_def_scaled / r_epi_scaled

    # phi is the ratio of the frequency of the epicycle to the frequency of the deferent.
    # The frequency of the nth harmonic is n times the fundamental frequency.
    phi = n_epi / n_def

    return R, phi

# Calculate the ordered pair (R, phi)
R, phi = get_R_phi()

# Print the final result in the specified format
# The values are integers, so we format them as such.
print(f"R = {R_def_scaled:.2f} / {r_epi_scaled:.2f} = {int(R)}")
print(f"phi = {n_epi} / {n_def} = {int(phi)}")
print(f"The ordered pair (R, phi) is: ({int(R)}, {int(phi)})")