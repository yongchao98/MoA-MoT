# --- Define Variables ---
# Let 'p' be the true per-exposure vaccine efficacy (VE_p) for an all-or-nothing vaccine.
# This means a proportion 'p' of the vaccinated population is fully protected.
p = 0.90 # Example: 90% efficacy, meaning 90% of recipients are fully immune.

# Let 'lambda_force' represent the force of infection (the hazard rate of infection
# for susceptible individuals). We can use a placeholder value, as it will cancel out.
lambda_force = 0.01 # Example: 1% chance of infection per unit of time for a susceptible person.

# --- Step-by-step Derivation ---
print("This script demonstrates whether 1 - Incidence Rate Ratio (IRR) correctly estimates")
print("the per-exposure vaccine efficacy (VE_p) for an all-or-nothing vaccine.\n")

print("Step 1: Define terms for an all-or-nothing vaccine.")
print(f"Let 'p' be the true per-exposure VE. This means a proportion 'p' of vaccinated individuals are fully immune.")
print(f"Let's assume a true VE, p = {p:.2f}")
# The remaining proportion of the vaccinated group is not protected.
proportion_susceptible = 1 - p
print(f"The remaining proportion, 1 - p = {proportion_susceptible:.2f}, remains fully susceptible.\n")

print("Step 2: Define Incidence Rates (IR).")
print("Let 'λ' be the force of infection (rate at which susceptibles get infected).")
# In the unvaccinated group, all individuals are susceptible.
# So, the incidence rate in the unvaccinated (IR_u) is equal to the force of infection.
ir_u = lambda_force
print(f"Incidence Rate in Unvaccinated (IR_u) = λ = {ir_u}")

# In the vaccinated group, only the proportion (1-p) is susceptible to infection.
# The number of cases is proportional to λ * (1-p).
# The person-time denominator for the rate calculation includes the entire group (both protected and unprotected).
# Therefore, the overall incidence rate for the vaccinated group is the force of infection scaled by the susceptible proportion.
ir_v = lambda_force * proportion_susceptible
print(f"Incidence Rate in Vaccinated (IR_v) = λ * (1 - p) = {lambda_force} * {proportion_susceptible:.2f} = {ir_v:.4f}\n")

print("Step 3: Calculate the Incidence Rate Ratio (IRR).")
# IRR is the ratio of the two incidence rates.
irr = ir_v / ir_u
print(f"IRR = IR_v / IR_u")
print(f"IRR = {ir_v:.4f} / {ir_u} = {irr:.2f}\n")

print("Step 4: Calculate Vaccine Efficacy from IRR.")
# The standard formula for VE from an IRR is 1 - IRR.
ve_from_irr = 1 - irr
print(f"Calculated VE = 1 - IRR")
# Outputting each number in the final equation:
print(f"Calculated VE = 1 - {irr:.2f} = {ve_from_irr:.2f}\n")

print("Step 5: Compare the calculated VE with the true VE.")
print(f"The calculated VE based on the IRR is {ve_from_irr:.2f}.")
print(f"The true per-exposure VE ('p') was defined as {p:.2f}.")
print("The calculated value is equal to the true value.\n")

print("--- Conclusion ---")
print("For an all-or-nothing vaccine, the formula 1 - IRR correctly estimates the per-exposure vaccine efficacy.")