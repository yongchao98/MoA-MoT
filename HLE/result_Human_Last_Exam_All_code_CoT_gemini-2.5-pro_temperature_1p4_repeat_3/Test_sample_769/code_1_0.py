def analyze_compliance():
    """
    Analyzes employers based on Ontario's employment laws for policies
    on disconnecting from work and electronic monitoring as of Jan 2, 2023.
    """
    print("--- Ontario Employment Law Compliance Analysis (as of Jan 2, 2023) ---\n")
    print("Applicable Rules:")
    print("1. An employer must have a written policy if they employ 25 or more employees on January 1st of a year.")
    print("2. For 2022, the 'Disconnecting from Work' policy was due by June 2, 2022.")
    print("3. For 2022, the 'Electronic Monitoring' policy was due by October 11, 2022.\n")

    employers = {
        'A': {'employees_jan_1_2022': 20, 'has_disconnect_policy': False, 'has_monitoring_policy': False},
        'B': {'employees_jan_1_2022': 23, 'has_disconnect_policy': True, 'has_monitoring_policy': True, 'distributes_to_new': False},
        'C': {'employees_jan_1_2022': 1000, 'has_disconnect_policy': True, 'has_monitoring_policy': True},
        'D': {'employees_jan_1_2022': 30, 'has_disconnect_policy': True, 'has_monitoring_policy': False},
        'E': {'employees_jan_1_2022': 22, 'has_disconnect_policy': False, 'has_monitoring_policy': False}
    }

    non_compliant_employer = None
    final_equation = ""

    for name, data in employers.items():
        print(f"--- Analyzing Employer {name} ---")
        
        count = data['employees_jan_1_2022']
        threshold = 25
        
        print(f"Number of employees on Jan 1, 2022: {count}")
        print(f"Compliance threshold: {threshold} employees")
        
        if count >= threshold:
            print(f"Evaluation: {count} >= {threshold}. The employer was required to have both policies.")
            
            # Check for missing policies
            if not data['has_disconnect_policy']:
                print("Result: NOT IN COMPLIANCE. Missing 'Disconnecting from Work' policy.")
                non_compliant_employer = name
                final_equation = f"The non-compliance is based on the calculation: {count} (employees) >= {threshold} (threshold), which triggers the policy requirement that was not met."

            elif not data['has_monitoring_policy']:
                print("Result: NOT IN COMPLIANCE. Missing 'Electronic Monitoring' policy.")
                non_compliant_employer = name
                final_equation = f"The non-compliance is based on the calculation: {count} (employees) >= {threshold} (threshold), which triggers the policy requirement that was not met."
            
            # Check for distribution issues (for Employer B ambiguity)
            elif 'distributes_to_new' in data and not data['distributes_to_new']:
                print("Result: NOT IN COMPLIANCE. Fails to distribute policy to new employees.")
                non_compliant_employer = name
                final_equation = f"The non-compliance is based on the calculation: {count} (employees) >= {threshold} (threshold), which triggers policy distribution requirements that were not met."
                
            else:
                print("Result: In Compliance. All requirements met.")
        else:
            print(f"Evaluation: {count} < {threshold}. The employer was NOT required to have these policies.")
            # Handle Employer B's ambiguity: if they were not required, their policies are voluntary.
            if name == 'B':
                print("Policies were likely voluntary, so failure to distribute may not be a violation of this specific law.")
            print("Result: In Compliance.")
            
        print("-" * 25 + "\n")

    # Final summary based on the most unambiguous case
    print("--- Final Conclusion ---")
    print("Employer D is unambiguously not in compliance.")
    print("They had 30 employees on Jan 1, 2022, which is over the threshold of 25, but failed to create a mandatory electronic monitoring policy.\n")
    
    # Final output with equation as requested
    print("Final Equation for Non-Compliant Employer (D):")
    non_compliant_count = employers['D']['employees_jan_1_2022']
    threshold = 25
    print(f"Number of Employees: {non_compliant_count}")
    print(f"Legal Threshold: {threshold}")
    print(f"The legal requirement was triggered because {non_compliant_count} >= {threshold}.")


analyze_compliance()
<<<D>>>