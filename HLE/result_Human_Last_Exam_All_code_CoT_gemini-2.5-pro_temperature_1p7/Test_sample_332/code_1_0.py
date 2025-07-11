# Step 1: Define initial parameters
lambda_ods1 = 500  # FIT
lambda_ods2 = 400  # FIT
lambda_ods3 = 700  # FIT
beta = 0.1
lambda_target = 100  # FIT for ASIL C
t_mission = 10000  # hours
fit_conversion = 1e-9  # Failures per hour for 1 FIT

# Step 2: Calculate independent failure rates
lambda_ods1_ind = (1 - beta) * lambda_ods1
lambda_ods2_ind = (1 - beta) * lambda_ods2
lambda_ods3_ind = (1 - beta) * lambda_ods3

print("--- Step 1: Independent and Common Cause Failure Rates ---")
print(f"Independent failure rate of ODS1 (λ_ODS1_ind): {lambda_ods1_ind:.2f} FIT")
print(f"Independent failure rate of ODS2 (λ_ODS2_ind): {lambda_ods2_ind:.2f} FIT")
print(f"Independent failure rate of ODS3 (λ_ODS3_ind): {lambda_ods3_ind:.2f} FIT")
print("")

# Step 3: Calculate the failure rate of the 2oo3 redundant part (from independent failures)
# Convert FIT to failures/hour for the formula
l1_h = lambda_ods1_ind * fit_conversion
l2_h = lambda_ods2_ind * fit_conversion
l3_h = lambda_ods3_ind * fit_conversion

# The average failure rate for a 2oo3 system is λ_avg ≈ (λ1*λ2 + λ1*λ3 + λ2*λ3) * t_mission
lambda_2oo3_h = (l1_h * l2_h + l1_h * l3_h + l2_h * l3_h) * t_mission

# Convert the result back to FIT
lambda_2oo3_fit = lambda_2oo3_h / fit_conversion

print("--- Step 2: 2-out-of-3 Redundant System Failure Rate ---")
print(f"Failure rate of the 2oo3 block (λ_2oo3_independent): {lambda_2oo3_fit:.2f} FIT")
print("")

# Step 4: Calculate the common cause failure rate (λ_ccf)
lambda_ods_avg = (lambda_ods1 + lambda_ods2 + lambda_ods3) / 3
lambda_ccf = beta * lambda_ods_avg

print("--- Step 3: Common Cause Failure (CCF) Rate ---")
print(f"Average failure rate of an ODS unit: {lambda_ods_avg:.2f} FIT")
print(f"CCF rate (λ_common_cause = β * λ_avg): {lambda_ccf:.2f} FIT")
print("")

# Step 5: Calculate the required voter failure rate
# λ_voter < λ_target - (λ_2oo3_fit + λ_ccf)
required_lambda_voter = lambda_target - (lambda_2oo3_fit + lambda_ccf)

print("--- Step 4: Final Calculation for Voter ---")
print("The total system failure rate is the sum of its parts:")
print(f"λ_system = λ_2oo3_independent + λ_common_cause + λ_voter")
print(f"To meet the target, λ_voter must be less than λ_target - (λ_2oo3_independent + λ_common_cause)")
print("\nFinal equation with calculated values:")
print(f"λ_voter < {lambda_target} - ({lambda_2oo3_fit:.2f} + {lambda_ccf:.2f})")
print(f"λ_voter < {lambda_target} - {(lambda_2oo3_fit + lambda_ccf):.2f}")
print(f"\nResult:")
print(f"λ_voter < {required_lambda_voter:.2f} FIT")
