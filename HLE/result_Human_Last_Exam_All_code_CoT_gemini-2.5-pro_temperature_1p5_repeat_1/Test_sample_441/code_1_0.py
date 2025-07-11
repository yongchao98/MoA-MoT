import math

# This script calculates the parameters (R, phi) for a deferent-epicycle model
# approximating an object's motion around a square at constant speed.

# R is the ratio of the deferent radius to the epicycle radius.
# phi is the ratio of the epicycle frequency to the deferent frequency.

# From Fourier analysis of the motion, we find that the deferent corresponds
# to the n=-3 harmonic and the epicycle to the n=1 harmonic. The radii are
# proportional to the magnitudes of the Fourier coefficients.
# Deferent radius |c_{-3}| is proportional to (3*pi+2)/9
# Epicycle radius |c_1| is proportional to (pi-2)

pi = math.pi

# Calculate the numerator and denominator for R
R_numerator = 3 * pi + 2
R_denominator = 9 * (pi - 2)

# Calculate R
R = R_numerator / R_denominator

# The deferent frequency is omega_d = -3 * omega_0
# The epicycle frequency is omega_e = 1 * omega_0
# phi is the ratio of these frequencies.
phi_numerator = 1
phi_denominator = -3
phi = phi_numerator / phi_denominator

# Output the results, showing the formulas as requested.
print("The parameter R is the ratio of the radii of the deferent and the epicycle.")
print(f"The equation for R is: R = (3 * pi + 2) / (9 * (pi - 2))")
print(f"Substituting pi = {pi}:")
print(f"R = ({R_numerator}) / ({R_denominator})")
print(f"R = {R}")
print("\nThe parameter phi is the ratio of the frequencies of the epicycle and the deferent.")
print(f"The equation for phi is: phi = (1 * omega_0) / (-3 * omega_0)")
print(f"phi = {phi}")
print(f"\nThe ordered pair is (R, phi) = ({R}, {phi})")
<<<({(3*math.pi+2)/(9*(math.pi-2))}, {-1/3})>>>