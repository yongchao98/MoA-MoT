import math

# --- Input Parameters ---
lambda_ods1_fit = 500  # FIT
lambda_ods2_fit = 400  # FIT
lambda_ods3_fit = 700  # FIT
beta = 0.10            # Common cause factor
system_target_fit = 100 # ASIL C target in FIT
t_mission_h = 10000    # Sample run time in hours

# --- Calculation ---

print("Step 1: Calculate the Common Cause Failure (CCF) rate (λ_ccf).")
# For non-identical components, we use the average failure rate.
lambda_avg_fit = (lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit) / 3
lambda_ccf_fit = beta * lambda_avg_fit
print(f"  Average ODS failure rate λ_avg = ({lambda_ods1_fit} + {lambda_ods2_fit} + {lambda_ods3_fit}) / 3 = {lambda_avg_fit:.2f} FIT")
print(f"  CCF rate λ_ccf = β * λ_avg = {beta} * {lambda_avg_fit:.2f} = {lambda_ccf_fit:.2f} FIT\n")

print("Step 2: Calculate the independent failure rate for each ODS (λ_ind).")
# The independent failure rate is the portion not subject to common cause failures.
lambda_1_ind_fit = (1 - beta) * lambda_ods1_fit
lambda_2_ind_fit = (1 - beta) * lambda_ods2_fit
lambda_3_ind_fit = (1 - beta) * lambda_ods3_fit
print(f"  λ_ind_1 = (1 - {beta}) * {lambda_ods1_fit} = {lambda_1_ind_fit:.0f} FIT")
print(f"  λ_ind_2 = (1 - {beta}) * {lambda_ods2_fit} = {lambda_2_ind_fit:.0f} FIT")
print(f"  λ_ind_3 = (1 - {beta}) * {lambda_ods3_fit} = {lambda_3_ind_fit:.0f} FIT\n")

print("Step 3: Calculate the failure rate from pairs of independent failures (λ_2oo3_ind).")
# This is approximated by λ_2oo3_ind ≈ Σ(λi*λj) * T_mission for i!=j
# The result needs to be converted from FIT^2 * hours to FIT.
# 1 FIT = 1e-9 / hour. So the conversion factor is 1e-9.
prod_12 = lambda_1_ind_fit * lambda_2_ind_fit
prod_13 = lambda_1_ind_fit * lambda_3_ind_fit
prod_23 = lambda_2_ind_fit * lambda_3_ind_fit
sum_of_products = prod_12 + prod_13 + prod_23
lambda_2oo3_ind_fit = sum_of_products * t_mission_h * 1e-9
print(f"  Sum of products of pairs = ({lambda_1_ind_fit:.0f}*{lambda_2_ind_fit:.0f}) + ({lambda_1_ind_fit:.0f}*{lambda_3_ind_fit:.0f}) + ({lambda_2_ind_fit:.0f}*{lambda_3_ind_fit:.0f}) = {sum_of_products} FIT^2")
print(f"  λ_2oo3_ind = {sum_of_products} * {t_mission_h}h * 10^-9 = {lambda_2oo3_ind_fit:.2f} FIT\n")

print("Step 4: Calculate the total failure rate of the 2oo3 ODS subsystem (λ_ODS_subsystem).")
lambda_ods_subsystem_fit = lambda_ccf_fit + lambda_2oo3_ind_fit
print(f"  λ_ODS_subsystem = λ_ccf + λ_2oo3_ind = {lambda_ccf_fit:.2f} FIT + {lambda_2oo3_ind_fit:.2f} FIT = {lambda_ods_subsystem_fit:.2f} FIT\n")

print("Step 5: Calculate the required failure rate for the voter (λ_voter).")
# The total system failure rate must be less than the ASIL C target.
# λ_ODS_subsystem + λ_voter < system_target_fit
lambda_voter_max_fit = system_target_fit - lambda_ods_subsystem_fit
print(f"  The equation is: λ_ODS_subsystem + λ_voter < {system_target_fit} FIT")
print(f"  Substituting the value: {lambda_ods_subsystem_fit:.2f} FIT + λ_voter < {system_target_fit} FIT")
print(f"  Therefore: λ_voter < {system_target_fit} FIT - {lambda_ods_subsystem_fit:.2f} FIT")

print("\n--- Final Answer ---")
print(f"λvoter < {lambda_voter_max_fit:.2f} FIT")

# Final answer in the required format
final_answer = f"{lambda_voter_max_fit:.2f}"