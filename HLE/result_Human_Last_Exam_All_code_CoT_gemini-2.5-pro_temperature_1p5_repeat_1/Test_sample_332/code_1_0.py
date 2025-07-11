# Given parameters
lambda_ods1 = 500  # FIT
lambda_ods2 = 400  # FIT
lambda_ods3 = 700  # FIT
beta = 0.10
system_target_lambda = 100  # FIT for ASIL C
t_mission = 10000  # hours
fit_conversion = 1e-9  # Failures per hour per FIT

# Step 1: Calculate the Common Cause Failure (CCF) rate
# For non-identical components, we use the average failure rate.
avg_lambda = (lambda_ods1 + lambda_ods2 + lambda_ods3) / 3
lambda_ccf = beta * avg_lambda

print("Step 1: Calculate Common Cause Failure (CCF) Rate")
print(f"Average ODS failure rate = ({lambda_ods1} + {lambda_ods2} + {lambda_ods3}) / 3 = {avg_lambda:.2f} FIT")
print(f"λ_ccf = β * Average λ = {beta} * {avg_lambda:.2f} = {lambda_ccf:.2f} FIT\n")

# Step 2: Calculate the failure rate of the 2-out-of-3 system from independent failures
# First, find the independent failure rate for each ODS
lambda_ods1_ind = lambda_ods1 * (1 - beta)
lambda_ods2_ind = lambda_ods2 * (1 - beta)
lambda_ods3_ind = lambda_ods3 * (1 - beta)

print("Step 2: Calculate Independent System Failure Rate (λ_2oo3_ind)")
print(f"λ_ods1_ind = {lambda_ods1} * (1 - {beta}) = {lambda_ods1_ind:.2f} FIT")
print(f"λ_ods2_ind = {lambda_ods2} * (1 - {beta}) = {lambda_ods2_ind:.2f} FIT")
print(f"λ_ods3_ind = {lambda_ods3} * (1 - {beta}) = {lambda_ods3_ind:.2f} FIT")

# Convert independent FIT rates to failures/hour for the formula
l1_ind_h = lambda_ods1_ind * fit_conversion
l2_ind_h = lambda_ods2_ind * fit_conversion
l3_ind_h = lambda_ods3_ind * fit_conversion

# Calculate the combinations of dual failures (products of rates)
# The formula for λ_2oo3 is (λ1λ2 + λ2λ3 + λ3λ1) * t
term12 = l1_ind_h * l2_ind_h
term23 = l2_ind_h * l3_ind_h
term13 = l1_ind_h * l3_ind_h
sum_of_products = term12 + term23 + term13

# Calculate the final rate in failures/hour
lambda_2oo3_ind_h = sum_of_products * t_mission

# Convert the result back to FIT
lambda_2oo3_ind = lambda_2oo3_ind_h / fit_conversion

print(f"\nλ_2oo3_ind ≈ (λ_ind1*λ_ind2 + λ_ind2*λ_ind3 + λ_ind3*λ_ind1) * t")
print(f"λ_2oo3_ind ≈ ({lambda_2oo3_ind:.2f} FIT)\n")


# Step 3: Calculate the required failure rate for the voter
# λ_system = λ_2oo3_ind + λ_ccf + λ_voter
# λ_voter <= λ_target - λ_2oo3_ind - λ_ccf
lambda_voter_max = system_target_lambda - lambda_2oo3_ind - lambda_ccf

print("Step 3: Calculate the Maximum Voter Failure Rate")
print(f"λ_voter ≤ λ_target - λ_ccf - λ_2oo3_ind")
print(f"λ_voter ≤ {system_target_lambda} - {lambda_ccf:.2f} - {lambda_2oo3_ind:.2f}")
print(f"λ_voter ≤ {lambda_voter_max:.2f} FIT")

# Final result in the required format
print("\n---")
print(f"Final Answer:")
print(f"λvoter < {lambda_voter_max:.2f} FIT")
