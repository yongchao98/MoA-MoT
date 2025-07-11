import math

# ---
# Plan Execution
# ---

# 1. Define known constants
b = -0.09  # Basquin exponent

# 2. State the necessary assumption
# The ratio of endurance limit to ultimate tensile strength (sigma_e / sigma_UTS) is not provided.
# A typical engineering value for some materials is around 0.4, which we will use for this estimation.
sigma_e_over_UTS = 0.4

print("Step 1: Determine the fatigue life at the reference stress level (σ_e).")
# From the Basquin relationship and the given anchor point (N=0.5 at sigma=sigma_UTS),
# we can derive the formula for cycles to failure (N_f) as: N_f = 0.5 * (σ_a / σ_UTS)^(1/b)
# We now calculate N_f1, which is the life at σ_a = σ_e.

one_over_b = 1.0 / b
N_f1 = 0.5 * (sigma_e_over_UTS)**(one_over_b)

print(f"Basquin exponent (b) = {b}")
print(f"Assumed σ_e / σ_UTS ratio = {sigma_e_over_UTS}")
print(f"Cycles to failure at σ_e (N_f1) = 0.5 * ({sigma_e_over_UTS})^(1/{b}) = {N_f1:.3f}")
print("-" * 50)

print("Step 2: Apply the Palmgren-Miner linear damage rule.")
# The total life N_life can be calculated by rearranging the Palmgren-Miner equation:
# N_life = N_f1 / [0.7 + 0.2 * (1.1)^(-1/b) + 0.1 * (1.2)^(-1/b)]
# Let's calculate each term in the denominator.

term1 = 0.70
term2 = 0.20 * (1.1)**(-one_over_b)
term3 = 0.10 * (1.2)**(-one_over_b)
damage_factor = term1 + term2 + term3

print("The fatigue life equation is: N_life = N_f1 / (Damage_at_1.0σ_e + Damage_at_1.1σ_e + Damage_at_1.2σ_e)")
print(f"Component 1 (from 70% life at 1.0*σ_e) = {term1:.4f}")
print(f"Component 2 (from 20% life at 1.1*σ_e) = {term2:.4f}")
print(f"Component 3 (from 10% life at 1.2*σ_e) = {term3:.4f}")
print(f"Total Damage Factor (denominator) = {damage_factor:.4f}")
print("-" * 50)

print("Step 3: Calculate the final estimated fatigue life.")
N_life = N_f1 / damage_factor
print(f"N_life = {N_f1:.3f} / {damage_factor:.4f}")
print(f"Estimated Fatigue Life (N_life) = {N_life:.3f} cycles")
print("-" * 50)

# 4. Format the final answer as requested
formatted_answer = f"{N_life:.3e}"
print(f"The final answer in scientific notation (3 decimal places) is: {formatted_answer}")

print(f"<<<{formatted_answer}>>>")