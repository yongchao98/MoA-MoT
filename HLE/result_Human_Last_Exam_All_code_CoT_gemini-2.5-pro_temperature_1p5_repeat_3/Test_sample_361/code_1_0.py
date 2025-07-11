import numpy as np

# This script calculates the critical chemical potential (mu_c) for scalar condensation
# in a 5D Einstein-Gauss-Bonnet background at T=0, for an operator of dimension Delta=3.

# 1. Define the given parameters
lambda_gb = 0.1  # The Gauss-Bonnet coupling
# The known value of the critical chemical potential for lambda_gb = 0 from literature.
mu_c_0 = 2.255

# 2. Explain the theoretical basis
# The critical chemical potential mu_c scales inversely with the effective AdS radius L_eff.
# The ratio of the new effective radius to the old one (L_eff/L) is given by:
# sqrt( (1 + sqrt(1 - 4 * lambda_gb)) / 2 )
# Therefore, mu_c(lambda_gb) = mu_c(0) / (L_eff/L)

# 3. Calculate the scaling factor for the effective radius
# This factor represents L_eff / L.
scaling_factor_L_eff = np.sqrt((1 + np.sqrt(1 - 4 * lambda_gb)) / 2)

# 4. Calculate the new critical chemical potential
mu_c_final = mu_c_0 / scaling_factor_L_eff

# 5. Print the results clearly, showing each number in the equation.
print("The calculation for the critical chemical potential is based on the formula:")
print("mu_c(lambda_gb) = mu_c(0) / sqrt((1 + sqrt(1 - 4 * lambda_gb)) / 2)\n")

print(f"Given values:")
print(f"mu_c(0) = {mu_c_0}")
print(f"lambda_gb = {lambda_gb}\n")

print(f"Step-by-step calculation:")
print(f"mu_c({lambda_gb}) = {mu_c_0} / sqrt((1 + sqrt(1 - 4 * {lambda_gb})) / 2)")
print(f"mu_c({lambda_gb}) = {mu_c_0} / sqrt((1 + sqrt({1 - 4 * lambda_gb})) / 2)")
print(f"mu_c({lambda_gb}) = {mu_c_0} / sqrt((1 + {np.sqrt(1 - 4 * lambda_gb)}) / 2)")
print(f"mu_c({lambda_gb}) = {mu_c_0} / {scaling_factor_L_eff}")
print(f"mu_c({lambda_gb}) = {mu_c_final}\n")

print(f"The final value of the critical chemical potential is:")
print(f"{mu_c_final}")
<<<2.3939223791241193>>>