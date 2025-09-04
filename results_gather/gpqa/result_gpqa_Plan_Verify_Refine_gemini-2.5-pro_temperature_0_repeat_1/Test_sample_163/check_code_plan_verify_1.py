# Given values for System 1
P1 = 2.0  # years
K1a = 10.0  # km/s
K1b = 5.0  # km/s

# Given values for System 2
P2 = 1.0  # year
K2a = 15.0  # km/s
K2b = 10.0  # km/s

# The total mass M is proportional to P * (Ka + Kb)^3.
# We want to find the ratio M1 / M2.

# Calculate the sum of the velocity amplitudes for each system
K1_sum = K1a + K1b
K2_sum = K2a + K2b

# Calculate the ratio of the masses
# M1 / M2 = (P1 * (K1_sum)^3) / (P2 * (K2_sum)^3)
mass_ratio = (P1 / P2) * (K1_sum / K2_sum)**3

print(f"The ratio of the mass of system_1 to system_2 is: {mass_ratio:.3f}")