import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

# --- Explanation and Plan ---
# For an all-or-nothing vaccine, a proportion 'p' of vaccinees are fully protected,
# and '(1-p)' are not protected at all. The true efficacy is 'p'.
# We will model the incidence rates to see what '1 - IRR' calculates.
#
# Plan:
# 1. Define the true vaccine efficacy 'p' (e.g., 80% or 0.8).
# 2. Assume an incidence rate for the unvaccinated group (lambda_U).
# 3. Calculate the incidence rate for the vaccinated group (lambda_V). Since only
#    the '(1-p)' proportion is at risk, lambda_V = (1-p) * lambda_U.
# 4. Calculate the Incidence Rate Ratio (IRR) = lambda_V / lambda_U.
# 5. Calculate the estimated VE = 1 - IRR.
# 6. Compare the estimated VE to the true 'p' and print the conclusion.

# --- Simulation ---

# Step 1: Define true 'all-or-nothing' vaccine efficacy (p)
true_ve = 0.80

# Step 2: Define a hypothetical incidence rate in the unvaccinated group (lambda_U)
# (e.g., 10 cases per 100 person-years)
lambda_unvaccinated = 0.10

print("--- Scenario Setup ---")
print(f"True 'All-or-Nothing' Vaccine Efficacy (p) = {true_ve:.2f}")
print(f"Incidence Rate in Unvaccinated (λ_U) = {lambda_unvaccinated:.2f}")
print("-" * 25)

# Step 3: Calculate the theoretical incidence rate in the vaccinated group (lambda_V)
lambda_vaccinated = (1 - true_ve) * lambda_unvaccinated

# Step 4: Calculate the Incidence Rate Ratio (IRR)
irr = lambda_vaccinated / lambda_unvaccinated

# Step 5: Calculate the vaccine efficacy estimate from the IRR
estimated_ve = 1 - irr

# --- Final Equation and Result ---
print("The final calculation is: Estimated VE = 1 - (λ_V / λ_U)")
print(f"Estimated VE = 1 - ({lambda_vaccinated:.4f} / {lambda_unvaccinated:.2f})")
print(f"Estimated VE = 1 - {irr:.2f}")
print(f"Estimated VE = {estimated_ve:.2f}")
print("-" * 25)

print("\n--- Conclusion ---")
print(f"The true vaccine efficacy was set to: {true_ve:.2f}")
print(f"The estimated vaccine efficacy (1 - IRR) is: {estimated_ve:.2f}")

if estimated_ve == true_ve:
    print("\nResult: 1 - IRR correctly estimates the per-exposure vaccine efficacy.")
else:
    print("\nResult: 1 - IRR does not correctly estimate the per-exposure vaccine efficacy.")

# Restore original stdout
sys.stdout = original_stdout
# Print the captured output to the console
print(buffer.getvalue())
<<<C>>>