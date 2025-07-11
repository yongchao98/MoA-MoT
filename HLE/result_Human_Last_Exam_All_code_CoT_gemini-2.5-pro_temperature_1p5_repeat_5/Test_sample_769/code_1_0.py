def check_ontario_employment_law_compliance():
    """
    Analyzes several employer scenarios to determine compliance with Ontario's
    'Disconnecting from Work' and 'Electronic Monitoring' policy laws
    as of January 2, 2023.
    """
    
    # Define the legal rule from Ontario's Employment Standards Act
    EMPLOYEE_THRESHOLD = 25
    
    # Structure data for each employer based on the question
    employers = {
        'A': { 'employees_on_jan_1_2022': 20, 'has_monitoring_policy': False },
        'B': { 'employees_on_jan_1_2022': 23, 'has_monitoring_policy': True }, # Not required, but has it
        'C': { 'employees_on_jan_1_2022': 1000, 'has_monitoring_policy': True },
        'D': { 'employees_on_jan_1_2022': 30, 'has_monitoring_policy': False },
        'E': { 'employees_on_jan_1_2022': 22, 'has_monitoring_policy': False }
    }
    
    print("--- Analysis of Ontario Employment Law Compliance ---")
    print(f"The law requires employers with {EMPLOYEE_THRESHOLD} or more employees on January 1st of a year to have written policies for:")
    print("1. Disconnecting from Work")
    print("2. Electronic Monitoring\n")
    print("The compliance check date is January 2, 2023, so we check the employee count from January 1, 2022.\n")
    
    non_compliant_employer = None
    
    for name, data in employers.items():
        count = data['employees_on_jan_1_2022']
        
        print(f"Analyzing Employer {name}...")
        print(f"  - Employees on Jan 1, 2022: {count}")
        
        if count >= EMPLOYEE_THRESHOLD:
            print(f"  - Result of check: {count} >= {EMPLOYEE_THRESHOLD}. The employer IS subject to the law.")
            # Scenario D is the only one in this bucket that lacks a required policy.
            if not data['has_monitoring_policy']:
                print("  - Compliance Status: NOT COMPLIANT. The employer was required to have an electronic monitoring policy but has not developed one.")
                non_compliant_employer = name
            else: # This applies to Employer C
                print("  - Compliance Status: COMPLIANT. The employer has the required policies in place.")
        else:
            print(f"  - Result of check: {count} < {EMPLOYEE_THRESHOLD}. The employer IS NOT subject to the law.")
            print("  - Compliance Status: COMPLIANT.")
            
        print("-" * 20)
        
    if non_compliant_employer:
        final_count = employers[non_compliant_employer]['employees_on_jan_1_2022']
        print(f"\nFinal Conclusion: Employer {non_compliant_employer} is not in compliance.")
        print("Final Equation:")
        print(f"  Number of Employees ({final_count}) >= Threshold ({EMPLOYEE_THRESHOLD}) -> Policy is Required")
        print(f"  Electronic Monitoring Policy Exists? -> No")
        print(f"  Outcome: Non-Compliant")


check_ontario_employment_law_compliance()
print("<<<D>>>")