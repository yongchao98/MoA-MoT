import math

# --- Step 1: Define the given constants and loading profile ---
# Basquin exponent
b = -0.09
# Assumed cycles for fatigue endurance limit (a standard engineering convention)
N_e = 1e6
# Cycle proportions for each stress level
c1 = 0.70  # at sigma_e
c2 = 0.20  # at 1.1 * sigma_e
c3 = 0.10  # at 1.2 * sigma_e

# --- Step 2: Calculate the terms for the Palmgren-Miner equation ---
# The fatigue life N_f is given by the formula:
# N_f = N_e / (c1 + c2 * (1.1**(-1/b)) + c3 * (1.2**(-1/b)))
# Let's calculate the damage factors r2 and r3, which are 1.1**(-1/b) and 1.2**(-1/b)

r2 = 1.1**(-1/b)
r3 = 1.2**(-1/b)

# Calculate the denominator of the equation, which represents the total damage per N_e cycles
denominator = c1 + c2 * r2 + c3 * r3

# --- Step 3: Calculate the final fatigue life ---
N_f = N_e / denominator

# --- Step 4: Print the components of the final equation and the result ---
print("This calculation estimates the fatigue life using the Palmgren-Miner rule.")
print("The governing equation is derived as: N_f = N_e / (c1 + c2*r2 + c3*r3)\n")

print(f"Assumed cycles for endurance limit, N_e = {N_e:.0e}")
print(f"Proportion of life at sigma_e, c1 = {c1:.2f}")
print(f"Proportion of life at 1.1*sigma_e, c2 = {c2:.2f}")
print(f"Proportion of life at 1.2*sigma_e, c3 = {c3:.2f}")
print("-" * 30)
print(f"Relative damage factor for 1.1*sigma_e, r2 = 1.1**(-1/b) = {r2:.3f}")
print(f"Relative damage factor for 1.2*sigma_e, r3 = 1.2**(-1/b) = {r3:.3f}")
print(f"Total damage denominator = {c1:.3f} + {c2:.3f}*{r2:.3f} + {c3:.3f}*{r3:.3f} = {denominator:.3f}")
print("-" * 30)
print(f"Calculated Fatigue Life, N_f = {N_e:.0e} / {denominator:.3f} = {N_f:,.0f} cycles")

# --- Step 5: Format the final answer as requested ---
exponent = math.floor(math.log10(N_f))
mantissa = N_f / (10**exponent)
final_answer_string = f"{mantissa:.3f}e{exponent}"

print(f"\nThe final estimated fatigue life in the format Aeb is: {final_answer_string}")

# Return final answer in specified format
print(f'<<<{final_answer_string}>>>')
