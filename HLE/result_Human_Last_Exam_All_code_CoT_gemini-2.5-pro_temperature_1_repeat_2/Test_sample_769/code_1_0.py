def analyze_ontario_employment_laws():
    """
    Analyzes several employer scenarios based on Ontario's employment laws
    regarding 'Disconnecting from Work' and 'Electronic Monitoring' policies
    as of January 2, 2023.
    """
    # --- Law Parameters ---
    EMPLOYEE_THRESHOLD = 25
    DETERMINATION_DATE = "January 1, 2022"
    ASSESSMENT_DATE = "January 2, 2023"
    MONITORING_POLICY_DEADLINE = "October 11, 2022"

    # --- Employer Data ---
    employers = {
        "A": {
            "count_jan1_2022": 20,
            "has_disconnect_policy": False,
            "has_monitoring_policy": False,
            "summary": "Had 20 employees on Jan 1, 2022. Did not have a disconnect policy."
        },
        "B": {
            "count_jan1_2022": 23,
            "has_disconnect_policy": True,
            "has_monitoring_policy": True,
            "summary": "Had 23 employees. Had both policies but refused to distribute to new staff."
        },
        "C": {
            "count_jan1_2022": 1000,
            "has_disconnect_policy": True,
            "has_monitoring_policy": True,
            "summary": "Had 1,000 employees. Had and properly distributed both policies."
        },
        "D": {
            "count_jan1_2022": 30,
            "has_disconnect_policy": True,
            "has_monitoring_policy": False,
            "summary": "Had 30 employees on Jan 1, 2022, but did not develop an electronic monitoring policy."
        },
        "E": {
            "count_jan1_2022": 22,
            "has_disconnect_policy": False,
            "has_monitoring_policy": False,
            "summary": "Had 22 employees on Jan 1, 2022. Did not have a disconnect policy."
        }
    }

    print(f"Analyzing employer compliance as of {ASSESSMENT_DATE}")
    print(f"The key rule: Employers with {EMPLOYEE_THRESHOLD} or more employees on {DETERMINATION_DATE} must have both a 'Disconnecting from Work' and an 'Electronic Monitoring' policy.\n")

    non_compliant_employer = None

    for name, data in employers.items():
        print(f"--- Analyzing Employer {name} ---")
        count = data['count_jan1_2022']
        print(f"Employee count on {DETERMINATION_DATE}: {count}")

        # Check if the law applies
        if count >= EMPLOYEE_THRESHOLD:
            print(f"Since {count} is >= {EMPLOYEE_THRESHOLD}, the employer was required to have both policies.")
            
            required_disconnect = data['has_disconnect_policy']
            required_monitoring = data['has_monitoring_policy']

            print(f"  - Has 'Disconnecting from Work' policy? {'Yes' if required_disconnect else 'No'}")
            print(f"  - Has 'Electronic Monitoring' policy? {'Yes' if required_monitoring else 'No'}")

            if required_disconnect and required_monitoring:
                print("Result: This employer appears to be in compliance.")
            else:
                print(f"Result: NOT in compliance. They failed to implement all required policies by the deadlines in 2022.")
                non_compliant_employer = name
        else:
            print(f"Since {count} is < {EMPLOYEE_THRESHOLD}, the employer was NOT required to have these policies.")
            print("Result: This employer is in compliance with these specific rules.")
        
        print("-" * (29 + len(name)))
        print("\n")

    if non_compliant_employer:
        print(f"Final Conclusion: Employer {non_compliant_employer} is not in compliance.")
        print(f"The reason is that they had {employers[non_compliant_employer]['count_jan1_2022']} employees on {DETERMINATION_DATE}, which is above the threshold of {EMPLOYEE_THRESHOLD}, but failed to develop a policy on electronic monitoring by the {MONITORING_POLICY_DEADLINE} deadline.")

    # The final answer format required by the user
    if non_compliant_employer:
        print(f"<<<{non_compliant_employer}>>>")

analyze_ontario_employment_laws()