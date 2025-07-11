def check_compliance_d():
    """
    Analyzes Employer D's compliance with Ontario employment laws
    regarding mandatory workplace policies.
    """
    # Scenario data for Employer D
    employer_id = 'D'
    employees_on_jan_1_2022 = 30
    has_disconnect_policy = True
    has_monitoring_policy = False
    
    # Legal rules in Ontario
    POLICY_THRESHOLD = 25
    
    print(f"Analysis for Employer {employer_id}")
    print("=" * 30)

    # Step 1: Check if the policy requirement applies based on employee count on the critical date.
    print(f"1. Checking employee count on January 1, 2022...")
    print(f"   Employer had {employees_on_jan_1_2022} employees.")
    print(f"   The legal threshold is {POLICY_THRESHOLD} employees.")
    
    is_policy_required = employees_on_jan_1_2022 >= POLICY_THRESHOLD
    
    print(f"\n   Is {employees_on_jan_1_2022} >= {POLICY_THRESHOLD}? {is_policy_required}")

    if not is_policy_required:
        print("\nResult: Employer is IN COMPLIANCE as they are below the threshold.")
        return

    print("\n   Conclusion: The employer is required to have BOTH policies.")
    
    # Step 2: Verify if the required policies exist.
    print("\n2. Verifying presence of required policies...")
    print(f"   Has 'Disconnecting from Work' policy? {has_disconnect_policy}")
    print(f"   Has 'Electronic Monitoring' policy? {has_monitoring_policy}")

    # Step 3: Final determination of compliance status.
    is_compliant = has_disconnect_policy and has_monitoring_policy
    
    print("\n" + "=" * 30)
    if is_compliant:
        print("Final Verdict: Employer is IN COMPLIANCE.")
    else:
        print("Final Verdict: Employer is NOT IN COMPLIANCE.")
        print("Reason: The employer met the employee threshold but failed to implement all required policies.")
        if not has_monitoring_policy:
            print("Specifically, the 'Electronic Monitoring' policy is missing.")
            
check_compliance_d()