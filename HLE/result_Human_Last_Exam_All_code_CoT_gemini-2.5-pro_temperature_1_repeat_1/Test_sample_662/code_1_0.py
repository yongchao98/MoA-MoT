import math

# --- Given Information ---
# Basic reproduction number
Ro = 3.0
# Percentage of vaccinated people who still get the infection
breakthrough_rate = 0.06

# --- Introduction ---
print("To prevent an epidemic, the effective reproduction number (Re) must be brought below 1.")
print("The critical vaccination coverage (p_c) is the proportion of the population that needs to be vaccinated to achieve this.")
print("The formula, which accounts for a vaccine that is not 100% effective, is: p_c = (1 - 1/Ro) / e")
print("where 'e' is the vaccine effectiveness.\n")

# --- Step 1: Calculate Vaccine Effectiveness (e) ---
# If 6% of vaccinated people get infected, the vaccine is effective at preventing infection in the other 94%.
vaccine_effectiveness = 1 - breakthrough_rate
print("Step 1: Calculate vaccine effectiveness (e) based on the breakthrough infection rate.")
print(f"Effectiveness (e) = 1 - Breakthrough Rate = 1 - {breakthrough_rate} = {vaccine_effectiveness}\n")

# --- Step 2: Calculate the Herd Immunity Threshold for a Perfect Vaccine ---
# This is the (1 - 1/Ro) part of the equation.
herd_immunity_threshold = 1 - (1 / Ro)
print("Step 2: Calculate the herd immunity threshold (the proportion of the population that needs to be immune).")
print(f"Herd Immunity Threshold = 1 - (1 / Ro) = 1 - (1 / {Ro:.1f}) = {herd_immunity_threshold:.4f}\n")

# --- Step 3: Calculate the Final Critical Vaccination Coverage (p_c) ---
# This is the herd immunity threshold adjusted for the vaccine's actual effectiveness.
critical_coverage = herd_immunity_threshold / vaccine_effectiveness
print("Step 3: Calculate the required vaccine coverage by adjusting for the vaccine's effectiveness.")
print("Final Equation:")
# The user wants to see each number in the final equation.
print(f"p_c = (1 - 1 / {Ro:.1f}) / {vaccine_effectiveness}")
print(f"p_c = {herd_immunity_threshold:.4f} / {vaccine_effectiveness}")
print(f"p_c = {critical_coverage:.4f}\n")

# --- Conclusion ---
print(f"Therefore, the theoretical threshold of vaccine coverage your state should reach is {critical_coverage:.1%}.")

# --- Final Answer ---
# The final answer is the percentage value, rounded to one decimal place.
final_answer = round(critical_coverage * 100, 1)
print(f'<<<__{final_answer}__>>>')