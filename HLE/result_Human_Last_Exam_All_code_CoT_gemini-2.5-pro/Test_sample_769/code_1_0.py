import sys

# Step 1: Define the legal requirements in Ontario.
# The compliance check is for Jan 2, 2023, so we look at the requirements based on the Jan 1, 2022 employee count.
EMPLOYEE_THRESHOLD = 25

# Step 2: Model the data for Employer D.
employer_d_employees_jan1_2022 = 30
employer_d_has_disconnect_policy = True
employer_d_has_monitoring_policy = False

# Step 3: Perform the compliance analysis.
print("Analyzing Compliance for Employer D as of January 2, 2023:")
print("-" * 55)

# Check if the legal requirement was triggered
was_requirement_triggered = employer_d_employees_jan1_2022 >= EMPLOYEE_THRESHOLD

print(f"1. Employee count on Jan 1, 2022: {employer_d_employees_jan1_2022}")
print(f"2. Legal threshold for policies: {EMPLOYEE_THRESHOLD}")
print(f"3. Was employer required to have policies? ({employer_d_employees_jan1_2022} >= {EMPLOYEE_THRESHOLD}) -> {was_requirement_triggered}")

if not was_requirement_triggered:
    print("\nResult: Employer was below the threshold and is compliant.")
    sys.exit()

print("\nEmployer was required to have TWO policies: 'Disconnecting from Work' and 'Electronic Monitoring'.")

# Check if the employer fulfilled all legal requirements
is_compliant = employer_d_has_disconnect_policy and employer_d_has_monitoring_policy

print("\n4. Checking fulfillment of policy requirements:")
print(f"   - Has 'Disconnecting from Work' Policy? -> {employer_d_has_disconnect_policy}")
print(f"   - Has 'Electronic Monitoring' Policy?  -> {employer_d_has_monitoring_policy}")

# Step 4: Conclude and present the final "equation" of the logic.
print("\nFinal Compliance Equation:")
print(f"({employer_d_has_disconnect_policy}) AND ({employer_d_has_monitoring_policy}) -> Compliance is {is_compliant}")

if not is_compliant:
    print("\nConclusion: Employer D is NOT in compliance with applicable employment laws.")
else:
    print("\nConclusion: Employer D is in compliance.")
