# Step 0: Define the given parameters
lambda_ods1 = 500  # FIT
lambda_ods2 = 400  # FIT
lambda_ods3 = 700  # FIT
beta = 0.10       # Common Cause Failure factor
lambda_target = 100 # FIT for ASIL C
T = 10000         # hours, mission time

# Conversion factor from FIT to failures per hour
FIT_to_per_hour = 1e-9

print("Step-by-step calculation for the required voter failure rate.\n")

# Step 1: Calculate the Common Cause Failure (CCF) rate (lambda_ccf)
# For non-identical components, we use the average failure rate of the group.
lambda_ods_avg = (lambda_ods1 + lambda_ods2 + lambda_ods3) / 3
lambda_ccf = beta * lambda_ods_avg
print(f"1. The average failure rate of the ODS components is {lambda_ods_avg:.2f} FIT.")
print(f"   The Common Cause Failure (CCF) rate is β * λ_avg = {beta} * {lambda_ods_avg:.2f} = {lambda_ccf:.2f} FIT.\n")

# Step 2: Calculate the independent failure rates of the ODS components
lambda_ods1_ind = (1 - beta) * lambda_ods1
lambda_ods2_ind = (1 - beta) * lambda_ods2
lambda_ods3_ind = (1 - beta) * lambda_ods3
print("2. The independent failure rates of the ODS components are:")
print(f"   λ'_ODS1 = (1 - {beta}) * {lambda_ods1} = {lambda_ods1_ind:.2f} FIT")
print(f"   λ'_ODS2 = (1 - {beta}) * {lambda_ods2} = {lambda_ods2_ind:.2f} FIT")
print(f"   λ'_ODS3 = (1 - {beta}) * {lambda_ods3} = {lambda_ods3_ind:.2f} FIT\n")

# Step 3: Calculate the failure rate of the 2-out-of-3 system from independent failures (lambda_2oo3_ind)
# Convert independent FIT rates to failures per hour for the formula
l1_ind_ph = lambda_ods1_ind * FIT_to_per_hour
l2_ind_ph = lambda_ods2_ind * FIT_to_per_hour
l3_ind_ph = lambda_ods3_ind * FIT_to_per_hour

# The formula for the average failure rate of a 2oo3 system over mission time T is:
# λ_2oo3 ≈ (λ'1*λ'2 + λ'1*λ'3 + λ'2*λ'3) * T
sum_of_products = (l1_ind_ph * l2_ind_ph) + \
                  (l1_ind_ph * l3_ind_ph) + \
                  (l2_ind_ph * l3_ind_ph)

lambda_2oo3_ind_ph = sum_of_products * T

# Convert the result back to FIT
lambda_2oo3_ind = lambda_2oo3_ind_ph / FIT_to_per_hour
print("3. The failure rate of the 2-out-of-3 ODS system due to independent failures is:")
print(f"   λ_2oo3_independent ≈ {lambda_2oo3_ind:.2f} FIT.\n")

# Step 4: Calculate the maximum allowed failure rate for the voter
# λ_voter < λ_target - λ_ccf - λ_2oo3_independent
lambda_voter_max = lambda_target - lambda_ccf - lambda_2oo3_ind

print("4. The final equation to find the required voter failure rate is:")
print("   λ_voter < λ_target - λ_ccf - λ_2oo3_independent")
print(f"   λ_voter < {lambda_target} FIT - {lambda_ccf:.2f} FIT - {lambda_2oo3_ind:.2f} FIT")
print(f"   λ_voter < {lambda_voter_max:.2f} FIT\n")

print("Final Answer:")
print(f"λvoter < {lambda_voter_max:.2f} FIT")