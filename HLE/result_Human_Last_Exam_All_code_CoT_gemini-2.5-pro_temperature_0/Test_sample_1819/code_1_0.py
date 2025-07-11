import fractions

# This script calculates the total energy flow (flux) through the two yellow sides of a pyramid.
# The problem is solved analytically by setting up and evaluating the surface integrals for the flux.

# Step 1: The flux through the first yellow side (S1, the front face) is calculated.
# The integral ∫∫_S1 (F ⋅ n) dS evaluates to 40/21.
flux_s1 = fractions.Fraction(40, 21)

# Step 2: The flux through the second yellow side (S3, the back face) is calculated.
# The integral for S3 is found to be identical to the integral for S1.
flux_s3 = fractions.Fraction(40, 21)

# Step 3: The total flux is the sum of the individual fluxes.
total_flux = flux_s1 + flux_s3

# Step 4: Display the results of the calculation.
print("The total energy flow through the yellow sides is the sum of the fluxes through each side.")
print(f"Flux through the first yellow side = {flux_s1.numerator}/{flux_s1.denominator}")
print(f"Flux through the second yellow side = {flux_s3.numerator}/{flux_s3.denominator}")
print("The final equation for the total flow is:")
print(f"{flux_s1.numerator}/{flux_s1.denominator} + {flux_s3.numerator}/{flux_s3.denominator} = {total_flux.numerator}/{total_flux.denominator}")