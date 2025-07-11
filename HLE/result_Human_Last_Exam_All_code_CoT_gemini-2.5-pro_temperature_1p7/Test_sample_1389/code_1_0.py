# 1. Define parameters for our model.
# Let's assume the vaccine provides perfect protection to 90% of recipients.
# This means the true per-exposure vaccine efficacy (VE_p) is 0.90.
per_exposure_VE = 0.90

# Let's assume an incidence rate in the unvaccinated group, e.g., 10 cases per 100 person-years.
IR_u = 0.10

print("Step 1: Define initial assumptions for an 'all-or-nothing' vaccine.")
print(f"True Per-Exposure Vaccine Efficacy (v): {per_exposure_VE}")
print(f"Incidence Rate in Unvaccinated Population (IR_u): {IR_u}\n")

# 2. Calculate the theoretical incidence rate in the vaccinated population (IR_v).
# In an all-or-nothing model, only the proportion (1 - v) of vaccinated people
# can get infected. Their rate of infection is the same as the unvaccinated.
# The other proportion (v) has a rate of 0.
# So, the overall rate is: IR_v = IR_u * (1 - v)
IR_v = IR_u * (1 - per_exposure_VE)

print("Step 2: Calculate the expected incidence rate in the vaccinated population (IR_v).")
print(f"IR_v = IR_u * (1 - v)")
print(f"IR_v = {IR_u} * (1 - {per_exposure_VE}) = {IR_v:.4f}\n")


# 3. Calculate the Incidence Rate Ratio (IRR) and the VE based on it.
IRR = IR_v / IR_u
VE_from_IRR = 1 - IRR

print("Step 3: Calculate the Vaccine Efficacy using the 1 - IRR formula.")
print("Equation: VE = 1 - (IR_v / IR_u)")
print(f"Final calculation: VE = 1 - ({IR_v:.4f} / {IR_u}) = {VE_from_IRR:.4f}\n")


# 4. Compare the result with the true per-exposure efficacy.
print("Step 4: Compare the calculated VE with the true per-exposure VE.")
print(f"True Per-Exposure VE (v) = {per_exposure_VE}")
print(f"Calculated VE from (1 - IRR) = {VE_from_IRR:.4f}\n")

print("Conclusion:")
print("The calculated value is identical to the true value.")
print("Therefore, 1 - IRR correctly estimates the per-exposure vaccine efficacy for an all-or-nothing vaccine.")

<<<C>>>