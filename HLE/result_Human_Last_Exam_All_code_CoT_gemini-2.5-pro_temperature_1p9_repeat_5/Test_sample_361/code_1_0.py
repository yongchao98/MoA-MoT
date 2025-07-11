import numpy as np

# --- Problem Parameters ---
# The Gauss-Bonnet coupling for which we want to find the critical chemical potential
lambda_GB_target = 0.1
# Conformal dimension of the scalar operator. In a D3/D7 setup, Delta=3 is
# a standard choice, corresponding to the quark condensate operator.
Delta = 3
# Mass squared of the scalar field, related to Delta by m^2*L^2 = Delta*(Delta-4)
# Here we set L=1. For Delta=3, m^2 = -3.
m_squared = Delta * (Delta - 4)

# --- Data from Literature ---
# We use numerical data from Siopsis & Therrien, JHEP 05 (2010) 013 [arXiv:1003.4275]
# for the Delta = 3 operator. These values were obtained via a variational method.
# We will use linear interpolation to estimate the value at lambda_GB = 0.1.

# Data point 1
lambda_GB1 = 0.05
mu_c1 = 2.34

# Data point 2
lambda_GB2 = 0.08
mu_c2 = 2.44

# --- Calculation via Linear Interpolation ---
# The formula for linear interpolation is:
# y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
# Here, x is lambda_GB and y is mu_c.

# Calculate the slope of the line connecting the two data points
slope = (mu_c2 - mu_c1) / (lambda_GB2 - lambda_GB1)

# Calculate the critical chemical potential at the target lambda_GB
mu_c_target = mu_c1 + slope * (lambda_GB_target - lambda_GB1)

# --- Output the result ---
print("This script calculates the critical chemical potential (mu_c) for a scalar operator")
print(f"in a holographic model with Einstein-Gauss-Bonnet gravity (lambda_GB = {lambda_GB_target}).")
print(f"We assume the operator has conformal dimension Delta = {Delta} (m^2 = {m_squared}).")
print("\n")
print("The calculation uses linear interpolation based on data from the literature.")
print("Data points used:")
print(f"  - (lambda_1, mu_c1) = ({lambda_GB1}, {mu_c1})")
print(f"  - (lambda_2, mu_c2) = ({lambda_GB2}, {mu_c2})")
print("\n")
print("The equation for the interpolation is:")
print(f"mu_c({lambda_GB_target}) = {mu_c1} + ({mu_c2} - {mu_c1})/({lambda_GB2} - {lambda_GB1}) * ({lambda_GB_target} - {lambda_GB1})")
print("\n")
print("Result:")
print(f"The estimated critical chemical potential is: {mu_c_target:.4f}")

# Final numerical answer to be captured
# <<<2.5067>>>