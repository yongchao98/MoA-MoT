import math

# Step 1: Define the indices 'n' for the dominant Fourier components.
# The deferent corresponds to the fundamental frequency, n=1.
# The most significant harmonic for the square path corresponds to n=-3.
n_deferent = 1
n_epicycle = -3

# Step 2: Calculate R, the ratio of the radii.
# The radius of each component is proportional to the amplitude of the Fourier coefficient |c_n|.
# For an object moving along the perimeter of a polygon, the coefficients |c_n| scale as 1/n^2.
# So, the ratio of radii R = |c_1| / |c_{-3}| can be calculated from the ratio of 1/n^2 terms.
R = (1/n_deferent**2) / (1/n_epicycle**2)
# R can be calculated more simply as:
R = (n_epicycle / n_deferent)**2

# Step 3: Calculate phi, the ratio of the frequencies.
# The deferent frequency is w_d = n_deferent * w0. For simplicity, we can let w0 = 1.
w_d = n_deferent
# The absolute frequency of the epicycle term is w_e_abs = n_epicycle * w0.
w_e_abs = n_epicycle

# In the historical deferent-epicycle model, the epicycle frequency w_ep is relative to
# the deferent arm. The relationship is w_e_abs = w_d + w_ep.
w_ep = w_e_abs - w_d

# phi is the ratio of the epicycle frequency (relative) to the deferent frequency.
phi = w_ep / w_d

# The ordered pair (R, phi)
ordered_pair = (R, phi)

# Output the results
print(f"The ratio of the deferent radius to the epicycle radius, R, is: {int(R)}")
print(f"The ratio of the epicycle frequency to the deferent frequency, phi, is: {int(phi)}")
print(f"The ordered pair (R, phi) is: {ordered_pair}")