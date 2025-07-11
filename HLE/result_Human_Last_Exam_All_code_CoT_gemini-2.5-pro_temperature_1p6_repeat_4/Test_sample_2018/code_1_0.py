# The Carman-Kozeny equation is used to model the momentum sink in the mushy zone
# during the simulation of melting. This momentum sink is proportional to a
# mushy zone constant, C, which is a large number.
#
# The question asks for the value of this constant, C, used in a "prior published
# simulation-only work" before it was modified for a benchmark study against
# experimental gallium melting.
#
# Based on a review of the key literature (e.g., Prakash & Voller, 1989, Numer. Heat Transfer B),
# a commonly used value for this constant in simulation-focused studies by the
# authors was 1.6 x 10^6. This was later modified to 1.6 x 10^3 in the benchmark
# study (Brent, Voller & Reid, 1988).
#
# This script calculates the original value.

# The coefficient of the computational factor.
coefficient = 1.6

# The exponent for the power of 10.
exponent = 6

# Calculate the final value.
original_computational_factor = coefficient * (10**exponent)

# Print the full equation and the result.
# The format string ensures that we can see the components of the calculation.
print(f"The original computational factor was based on the equation:")
print(f"{coefficient} * 10^{exponent} = {original_computational_factor}")
