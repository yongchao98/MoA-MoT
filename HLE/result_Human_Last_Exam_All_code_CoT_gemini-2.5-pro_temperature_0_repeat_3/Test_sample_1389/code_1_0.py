import sys

# Plan:
# 1. Define the parameters for an "all-or-nothing" vaccine model.
# 2. Calculate the incidence rate in the unvaccinated group (IR_u).
# 3. Calculate the incidence rate in the vaccinated group (IR_v).
# 4. Compute the Incidence Rate Ratio (IRR).
# 5. Calculate the estimated vaccine efficacy (VE) as 1 - IRR.
# 6. Compare the estimated VE with the true VE to determine if it's an overestimate, underestimate, or correct estimate.

# Step 1: Define parameters
# True vaccine efficacy (VE_S): The proportion of vaccinated people who are fully protected.
# Let's assume 90% efficacy.
VE_S = 0.90

# Force of infection (lambda): The rate at which susceptible individuals get infected.
# We can use an arbitrary value, as it will cancel out in the ratio.
# Let's say 5 cases per 100 person-years.
lambda_force_of_infection = 0.05

print("--- Model Parameters ---")
print(f"True Vaccine Efficacy (VE_S, proportion protected): {VE_S:.2f}")
print(f"Force of Infection (lambda): {lambda_force_of_infection}\n")

# Step 2: Calculate Incidence Rate in Unvaccinated (IR_u)
# In the unvaccinated group, everyone is assumed to be susceptible.
# So, the incidence rate is equal to the force of infection.
IR_u = lambda_force_of_infection
print("--- Incidence Rate Calculation ---")
print(f"Incidence Rate in Unvaccinated (IR_u) = lambda = {IR_u}")

# Step 3: Calculate Incidence Rate in Vaccinated (IR_v)
# In the vaccinated group, only the proportion (1 - VE_S) is susceptible.
# The incidence rate is the force of infection acting on this susceptible proportion.
proportion_susceptible_in_vax_group = 1 - VE_S
IR_v = lambda_force_of_infection * proportion_susceptible_in_vax_group
print(f"Incidence Rate in Vaccinated (IR_v) = lambda * (1 - VE_S) = {IR_v:.4f}\n")

# Step 4 & 5: Calculate Estimated VE from IRR
print("--- Efficacy Calculation from IRR ---")
# The Incidence Rate Ratio (IRR) is IR_v / IR_u
IRR = IR_v / IR_u

# The estimated vaccine efficacy is 1 - IRR.
VE_estimated = 1 - IRR

# Step 6: Show the final equation with numbers and conclude
print("The final equation for estimated VE is: VE_est = 1 - IRR")
print("Substituting the values:")
print(f"VE_est = 1 - (IR_v / IR_u)")
# Using sys.stdout.write to prevent the extra space from print's sep parameter
sys.stdout.write("VE_est = 1 - (")
sys.stdout.flush()
print(f"{IR_v:.4f} / {IR_u})")
print(f"VE_est = 1 - {IRR:.2f}")
print(f"VE_est = {VE_estimated:.2f}\n")

print("--- Conclusion ---")
print(f"The true vaccine efficacy (VE_S) was: {VE_S:.2f}")
print(f"The estimated vaccine efficacy (1 - IRR) is: {VE_estimated:.2f}")

if abs(VE_estimated - VE_S) < 1e-9:
    conclusion = "correctly estimate"
elif VE_estimated > VE_S:
    conclusion = "overestimate"
else:
    conclusion = "underestimate"

print(f"\nResult: For an all-or-nothing vaccine, 1 - IRR will {conclusion} the per-exposure vaccine efficacy.")
