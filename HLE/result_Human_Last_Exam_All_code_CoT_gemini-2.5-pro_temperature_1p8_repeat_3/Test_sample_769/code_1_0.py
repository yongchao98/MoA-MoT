def check_compliance_as_of_jan_2_2023():
    """
    Analyzes several employer scenarios to determine compliance with
    Ontario's 'disconnecting from work' and 'electronic monitoring' policy laws
    as of January 2, 2023.
    """
    print("Analyzing Ontario Employment Law Compliance as of January 2, 2023\n")
    print("Rule: An employer must have a 'disconnecting from work' policy AND an 'electronic monitoring' policy if they had 25 or more employees on January 1, 2022.\n")

    employers = {
        'A': {
            'description': 'Employed 20 employees on Jan 1, 2022. Did not develop a disconnect policy.',
            'count_jan1_2022': 20,
            'has_disconnect_policy': False,
            'has_monitoring_policy': False,
        },
        'B': {
            'description': 'Has 23 employees and both policies, but refuses to distribute them to new employees. The count on Jan 1, 2022 is assumed to be 23 as no other information is given.',
            'count_jan1_2022': 23,
            'has_disconnect_policy': True,
            'has_monitoring_policy': True,
        },
        'C': {
            'description': 'Had 1000 employees on Jan 1, 2022. Has and distributes both policies.',
            'count_jan1_2022': 1000,
            'has_disconnect_policy': True,
            'has_monitoring_policy': True,
        },
        'D': {
            'description': 'Employed 30 on Jan 1, 2022. Fired 6 employees in March 2022. Has a disconnect policy but not an electronic monitoring policy.',
            'count_jan1_2022': 30,
            'has_disconnect_policy': True,
            'has_monitoring_policy': False,
        },
        'E': {
            'description': 'Had 22 employees on Jan 1, 2022. Did not develop a disconnect policy.',
            'count_jan1_2022': 22,
            'has_disconnect_policy': False,
            'has_monitoring_policy': False,
        }
    }
    
    non_compliant_employer = None

    for name, data in employers.items():
        print(f"--- Checking Employer {name} ---")
        count = data['count_jan1_2022']
        
        print(f"Number of employees on January 1, 2022: {count}")
        
        # Determine if policies are required
        requires_policies = count >= 25
        
        if requires_policies:
            print(f"Result: With {count} employees, the employer was REQUIRED to have both policies.")
            # Check for compliance
            if data['has_disconnect_policy'] and data['has_monitoring_policy']:
                print("Status: Employer has both required policies. Compliant.")
            else:
                print("Status: Employer is MISSING one or more required policies. NOT IN COMPLIANCE.")
                if not data['has_disconnect_policy']:
                    print("Reason: Missing 'disconnecting from work' policy.")
                if not data['has_monitoring_policy']:
                    print("Reason: Missing 'electronic monitoring' policy.")
                non_compliant_employer = name
        else:
            print(f"Result: With {count} employees, the employer was NOT REQUIRED to have the policies.")
            print("Status: Compliant.")
        
        print("-" * (len(name) + 22) + "\n")

    if non_compliant_employer:
        print(f"Final Conclusion: Employer {non_compliant_employer} is not in compliance with applicable employment laws.")
    else:
        print("Final Conclusion: All employers evaluated appear to be in compliance based on the data.")

check_compliance_as_of_jan_2_2023()
<<<D>>>