def analyze_employer_compliance():
    """
    Analyzes employers for compliance with Ontario's employment laws
    regarding 'Disconnecting from Work' and 'Electronic Monitoring' policies
    as of January 2, 2023.
    """

    EMPLOYEE_THRESHOLD = 25

    # Data for each employer based on the question
    employers = {
        'A': {'count_jan_1_2022': 20, 'has_disconnect_policy': False, 'has_monitoring_policy': False, 'desc': '20 employees, no disconnect policy.'},
        'B': {'count_jan_1_2022': 23, 'has_disconnect_policy': True, 'has_monitoring_policy': True, 'desc': '23 employees, has policies but has issues with distribution.'},
        'C': {'count_jan_1_2022': 1000, 'has_disconnect_policy': True, 'has_monitoring_policy': True, 'desc': '1000+ employees, has and distributes both policies.'},
        'D': {'count_jan_1_2022': 30, 'has_disconnect_policy': True, 'has_monitoring_policy': False, 'desc': '30 employees, missing electronic monitoring policy.'},
        'E': {'count_jan_1_2022': 22, 'has_disconnect_policy': False, 'has_monitoring_policy': False, 'desc': '22 employees, no disconnect policy.'}
    }

    print("Analyzing Employer Compliance as of January 2, 2023")
    print(f"Legal Rule: Employers with {EMPLOYEE_THRESHOLD} or more employees on January 1, 2022, must have required policies.\n")

    non_compliant_employer = None

    for name, data in employers.items():
        count = data['count_jan_1_2022']
        
        print(f"--- Analysis for Employer {name} ---")
        print(f"Details: {data['desc']}")
        
        is_above_threshold = (count >= EMPLOYEE_THRESHOLD)
        # This is the "equation" part of the analysis
        print(f"Equation: Number of employees on Jan 1, 2022 ({count}) >= Threshold ({EMPLOYEE_THRESHOLD})? Result: {is_above_threshold}")

        if is_above_threshold:
            print("Requirement: Must have BOTH a 'Disconnecting from Work' and an 'Electronic Monitoring' policy.")
            
            # Check for disconnect policy
            req1 = data['has_disconnect_policy']
            print(f"Equation: 'Disconnect from Work' policy exists? Result: {req1}")

            # Check for monitoring policy
            req2 = data['has_monitoring_policy']
            print(f"Equation: 'Electronic Monitoring' policy exists? Result: {req2}")

            is_compliant = req1 and req2
            
            if not is_compliant:
                print("Conclusion: NOT IN COMPLIANCE\n")
                non_compliant_employer = name
            else:
                # Employer C has both required policies.
                print("Conclusion: COMPLIANT\n")
        else:
            print("Requirement: Not statutorily required to have the policies as the employee count is below the threshold.")
            # Employer B is under the threshold, so the statutory requirements to create or distribute policies do not apply.
            print("Conclusion: COMPLIANT\n")

    print("--- Final Result ---")
    if non_compliant_employer:
        print(f"The employer not in compliance is Employer {non_compliant_employer} because they met the employee threshold but failed to create a legally required policy.")
    else:
        print("All employers appear to be in compliance based on the analysis.")

if __name__ == '__main__':
    analyze_employer_compliance()