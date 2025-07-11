import math

# Based on the physical analysis, we determine the exponents n_1 through n_6.

# Step 1: Determine n1, n2, n3 from the zero-frequency limit.
# In the limit omega -> 0, the magnetic field noise power spectral density S_B is proportional
# to the source strength (proportional to conductivity sigma and temperature T) and a geometric factor.
# For a thin sheet of thickness t at distance z, S_B ~ sigma * T * t / z^4.
# The problem gives S_B ~ sigma^{n1} * T^{n2} * z^{n3}.
# By comparing the two expressions, we find:
n1 = 1
n2 = 1
n3 = -4

# Step 2: Determine n4 for the low-frequency regime (omega << 1/(sigma*z*t)).
# In this regime, inductive screening is negligible, and the noise spectrum is white.
# This means the power spectral density is independent of frequency.
# S_B(omega) ~ omega^{n4}, so n4 must be 0 for a frequency-independent spectrum.
n4 = 0

# Step 3: Determine n5 for the intermediate-frequency regime (1/(sigma*z*t) << omega << 1/(sigma*t^2)).
# In this regime, the slab is inductively screened. The power spectrum is suppressed by a factor
# of approximately 1 / (omega*sigma*z*t)^2.
# So, S_B(omega) ~ (sigma*T*t/z^4) / (omega*sigma*z*t)^2, which simplifies to S_B ~ T / (omega^2 * sigma * t * z^6).
# The scaling with frequency is omega^{-2}.
n5 = -2

# Step 4: Determine n6 for the high-frequency regime (omega >> 1/(sigma*t^2)).
# Here, the skin depth delta = sqrt(2/(omega*sigma)) is less than t.
# We must replace t with delta in our formulas for source strength and screening.
# Source strength ~ sigma*T*delta ~ T*sqrt(sigma/omega).
# Screening factor ~ 1 / (omega*sigma*z*delta)^2 ~ 1 / (z^2 * omega * sigma).
# S_B(omega) ~ (T*sqrt(sigma/omega)) / (z^2*omega*sigma) ~ T / (omega^{3/2} * sqrt(sigma)).
# The scaling with frequency is omega^{-3/2}.
n6 = -1.5

# Step 5: Calculate the final sum: sum_{k=1 to 6} k * n_k.
# We compute 1*n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6
total_sum = 1*n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6

# Print the final equation with the determined values and the final result.
print(f"The calculation is:")
print(f"1*n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6")
print(f"= 1*({n1}) + 2*({n2}) + 3*({n3}) + 4*({n4}) + 5*({n5}) + 6*({n6})")
print(f"= {1*n1} + {2*n2} + ({3*n3}) + {4*n4} + ({5*n5}) + ({6*n6})")
print(f"= {total_sum}")
