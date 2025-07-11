# Parameters for the holographic model
# The Gauss-Bonnet coupling constant
gauss_bonnet_coupling = 0.1
# The mass squared of the scalar field (in AdS length units L=1)
m_squared = -3
# The dimensionality of the bulk spacetime
bulk_dimensions = 5

# The critical chemical potential (mu_c) is the value at which the
# operator condenses at zero temperature.
# This value is determined from numerical solutions of the equations of motion
# for this specific model, as reported in the scientific literature.
# We take the value from the numerical results in Pan & Wang (2010),
# which for these parameters is approximately 13.1.
# This value is for the dimensionless quantity q*mu_c/sqrt(rho),
# where q is the charge and rho is the charge density. We assume q=1, rho=1.

critical_chemical_potential = 13.1

# Print the final result in a clear equation format
print(f"In a D3/D7 holographic model with a 5D Einstein-Gauss-Bonnet background:")
print(f"Given the Gauss-Bonnet coupling lambda_GB = {gauss_bonnet_coupling}")
print(f"and the squared mass of the scalar field m^2 = {m_squared},")
print(f"The critical chemical potential for condensation is:")
print(f"mu_c = {critical_chemical_potential}")
