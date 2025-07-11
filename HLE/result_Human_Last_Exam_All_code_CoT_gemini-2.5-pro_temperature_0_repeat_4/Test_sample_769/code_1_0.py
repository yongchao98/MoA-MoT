def analyze_employer_d_compliance():
    """
    Analyzes the compliance status of Employer D based on Ontario employment law.
    """
    # --- Legal and Scenario-Specific Variables ---
    EMPLOYEE_THRESHOLD = 25
    COUNT_DATE = "January 1, 2022"
    MONITORING_POLICY_DEADLINE = "October 11, 2022"

    employees_on_jan_1 = 30
    employees_fired = 6
    has_monitoring_policy = False

    print("--- Analyzing Compliance for Employer D ---")

    # Step 1: Check if the employer met the employee threshold on the key date.
    print(f"\n1. On {COUNT_DATE}, the law required employers with {EMPLOYEE_THRESHOLD} or more employees to create certain policies.")
    print(f"   Employer D had {employees_on_jan_1} employees on this date.")

    # Step 2: Determine the legal obligation based on the threshold.
    if employees_on_jan_1 >= EMPLOYEE_THRESHOLD:
        print(f"\n2. Since {employees_on_jan_1} is greater than or equal to {EMPLOYEE_THRESHOLD}, the employer was legally obligated.")
        print(f"   A key obligation was to have a written electronic monitoring policy by {MONITORING_POLICY_DEADLINE}.")
        # Note on headcount change
        print(f"   The subsequent termination of {employees_fired} employees does not change this obligation for the year 2022.")
    else:
        # This case does not apply to Employer D
        print(f"\n2. Since {employees_on_jan_1} is less than {EMPLOYEE_THRESHOLD}, the employer was not legally obligated.")
        return

    # Step 3: Verify if the employer fulfilled their obligation.
    print(f"\n3. Did the employer have the required policy? The answer is: {has_monitoring_policy}.")

    # Step 4: Final Conclusion
    print("\n--- Final Conclusion ---")
    if not has_monitoring_policy:
        print(f"Employer D is NOT in compliance with applicable employment laws.")
        print(f"The final equation for the requirement is: Employee count ({employees_on_jan_1}) >= Threshold ({EMPLOYEE_THRESHOLD}).")
        print("Because this condition was met, but the policy was not created, the employer is non-compliant.")
    else:
        print("Employer D is in compliance.")

# Run the analysis
analyze_employer_d_compliance()