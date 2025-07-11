# The Lindhard polarization function in the specified limit is given by:
# Chi(0,0) = C * (n/E_F), where n is the electron density and E_F is the Fermi energy.
# Our goal is to find the numerical value of the prefactor, C.

# From the first-principles derivation, the prefactor C is equal to -3/2.
# The numbers in the final equation for this prefactor are the numerator and denominator.
numerator = 3
denominator = 2

# Perform the calculation to find the numerical value.
result = -numerator / denominator

# Print the final equation with each number explicitly shown, along with the result.
print(f"The numerical prefactor is calculated from the expression -({numerator}/{denominator}).")
print(f"Final value: {result}")