import numpy as np

# The problem asks for the size of the set Omega of initial conditions (a(0), b(0))
# within the domain [-10, 1] x [10, 20] that lead to a specific type of blow-up.

# Based on the analysis of the ODE system, the blow-up condition
# (a(t) -> +inf and b(t) -> -inf) occurs if and only if a(0) > 0.

# The domain for the initial conditions is:
a0_domain = [-10, 1]
b0_domain = [10, 20]

# The subset Omega corresponds to initial conditions where a(0) > 0.
# So, the interval for a(0) for points in Omega is (0, 1].
# The interval for b(0) is the full range [10, 20].
a0_omega = (0, 1]
b0_omega = [10, 20]

# We calculate the measure (area) of Omega by finding the lengths of these intervals.
a0_omega_length = a0_omega[1] - a0_omega[0]
b0_omega_length = b0_omega[1] - b0_omega[0]

# The measure m(Omega) is the product of the lengths.
m_omega = a0_omega_length * b0_omega_length

print("Step-by-step calculation of the measure of Omega:")
print(f"1. The interval for a(0) in Omega is ({a0_omega[0]}, {a0_omega[1]}].")
print(f"   The length of this interval is {a0_omega[1]} - {a0_omega[0]} = {a0_omega_length}")
print(f"2. The interval for b(0) in Omega is [{b0_omega[0]}, {b0_omega[1]}].")
print(f"   The length of this interval is {b0_omega[1]} - {b0_omega[0]} = {b0_omega_length}")
print("3. The measure of Omega is the product of these lengths (area of the region).")
print(f"   m(Omega) = {a0_omega_length} * {b0_omega_length}")
print(f"Final result: m(Omega) = {m_omega}")