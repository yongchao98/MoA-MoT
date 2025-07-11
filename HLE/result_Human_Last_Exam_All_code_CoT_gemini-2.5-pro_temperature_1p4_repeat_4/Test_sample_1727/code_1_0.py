import numpy as np

# --- Introduction ---
# The problem asks for the stable amplitude factor chi, which relates the 
# vertical surface displacement (topography) to the basal shear stress driving the flow.
# The factor is derived by solving the Stokes equations for a viscous fluid layer.
# The final expression for chi depends on the wavenumber 'k' and the layer thickness 'h'.

# --- Parameters ---
# We will use plausible geophysical values to illustrate the calculation.
# Let's assume the harmonic shear stress has a wavelength of 100 km.
lambda_val = 100e3  # Wavelength in meters (100 km)
# Let's assume the tectonic plate has a thickness of 30 km.
h = 30e3         # Plate thickness in meters (30 km)

# --- Calculation ---
# 1. Calculate the wavenumber k from the wavelength lambda.
k = 2 * np.pi / lambda_val

# 2. Calculate the non-dimensional parameter k*h.
kh = k * h

# 3. Calculate the stable amplitude factor chi using the derived formula.
# The formula is: chi = sinh(k*h) + k*h * cosh(k*h)
chi = np.sinh(kh) + kh * np.cosh(kh)

# --- Output ---
print("The stable amplitude factor chi relates the surface topography e_s to the basal stress S_0 by:")
print("e_s = chi * (S_0 / (Delta_rho * g))")
print("\nThe derived formula for the factor chi is:")
print("chi = sinh(k*h) + k*h * cosh(k*h)")
print("\nFor a numerical example with plausible values:")
print(f"Wavelength (lambda) = {lambda_val/1000} km")
print(f"Plate thickness (h) = {h/1000} km")
print(f"Wavenumber (k) = {k:.6f} m^-1")

print("\nThe final equation with each number substituted is:")
# Using f-string to format the output string with calculated values
print(f"chi = sinh({kh:.4f}) + {kh:.4f} * cosh({kh:.4f})")
print(f"chi = {np.sinh(kh):.4f} + {kh:.4f} * {np.cosh(kh):.4f}")
print(f"chi = {np.sinh(kh):.4f} + {kh * np.cosh(kh):.4f}")
print(f"chi = {chi:.4f}")
print("\nTherefore, the stable amplitude factor is approximately {:.2f} for these parameters.".format(chi))
<<<9.55>>>