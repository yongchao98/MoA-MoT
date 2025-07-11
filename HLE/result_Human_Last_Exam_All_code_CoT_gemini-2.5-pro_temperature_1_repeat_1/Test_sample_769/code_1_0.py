import sys
from datetime import date

# Define the rules based on Ontario's Employment Standards Act (ESA)
EMPLOYEE_THRESHOLD = 25
POLICY_UPDATE_DISTRIBUTION_DEADLINE_DAYS = 30
ASSESSMENT_DATE = date(2023, 1, 2)
ELECTRONIC_MONITORING_POLICY_DEADLINE_2022 = date(2022, 10, 11)

def analyze_compliance(employer):
    """
    Analyzes an employer's compliance with Ontario's disconnecting from work
    and electronic monitoring policy laws as of January 2, 2023.
    """
    print(f"--- Analyzing Employer {employer['name']} ---")
    is_compliant = True
    reasons = []

    # Step 1: Check requirements based on Jan 1, 2022 employee count
    count_2022 = employer.get('count_jan_1_2022')
    print(f"Employee count on Jan 1, 2022: {count_2022}")

    if count_2022 >= EMPLOYEE_THRESHOLD:
        print(f"As the count ({count_2022}) was 25 or more, the employer was required to have:")
        print("  1. A 'Disconnecting from Work' policy by March 1, 2022.")
        print(f"  2. An 'Electronic Monitoring' policy by October 11, 2022.")

        # Check for Disconnecting from Work Policy
        if not employer.get('has_disconnect_policy'):
            is_compliant = False
            reasons.append("Failed to develop a 'Disconnecting from Work' policy required in 2022.")
        
        # Check for Electronic Monitoring Policy
        if not employer.get('has_monitor_policy'):
            is_compliant = False
            reasons.append("Failed to develop an 'Electronic Monitoring' policy by the Oct 11, 2022 deadline.")

        # Check for policy update and distribution compliance
        if employer.get('updates_monthly') and employer.get('distributes_quarterly'):
            # A quarterly distribution (~90 days) cannot meet a 30-day deadline for monthly updates.
            is_compliant = False
            reasons.append("Policy updates are monthly, but distribution is quarterly, which violates the 30-day distribution rule for changes.")

    else:
        print(f"As the count ({count_2022}) was less than 25, the employer was not required to have these policies in 2022.")

    # Step 2: Check requirements based on Jan 1, 2023 employee count
    # Note: For the assessment date of Jan 2, 2023, no one is yet non-compliant with 2023 rules,
    # as the deadline for 2023 is March 1, 2023. We only check for lingering non-compliance from 2022.
    count_2023 = employer.get('count_jan_1_2023')
    print(f"Employee count on Jan 1, 2023: {count_2023}")
    if count_2023 < EMPLOYEE_THRESHOLD and is_compliant:
         print(f"The employer is also not required to create new policies for 2023 (deadline March 1, 2023).")
    
    # Step 3: Handle specific cases like distribution refusal
    if employer.get('refuses_distribution'):
        # If the employer wasn't required to have the policy, refusing to distribute it isn't a violation
        # of these specific ESA sections.
        if count_2023 < EMPLOYEE_THRESHOLD:
            print("Although refusing to distribute policies is poor practice, it is not a violation of these specific laws as the employer is below the 25-employee threshold.")
        else: # This case isn't in the options, but for completeness:
            is_compliant = False
            reasons.append("Refusing to distribute required policies to employees is a violation.")
    
    # Step 4: Final conclusion for the employer
    print("\nConclusion:")
    if is_compliant:
        print(f"Employer {employer['name']} is IN COMPLIANCE with these specific laws as of Jan 2, 2023.")
    else:
        print(f"Employer {employer['name']} is NOT IN COMPLIANCE as of Jan 2, 2023.")
        for reason in reasons:
            print(f"  - Reason: {reason}")
    print("-" * 25 + "\n")
    return is_compliant, employer['name']


def solve_task():
    """
    Sets up the data for each employer and runs the analysis to find the non-compliant one.
    """
    employers = [
        {
            'name': 'A',
            'count_jan_1_2022': 20,
            'count_jan_1_2023': 20, # Hires start Jan 10, not counted on Jan 1
            'has_disconnect_policy': False,
            'has_monitor_policy': False,
        },
        {
            'name': 'B',
            'count_jan_1_2022': 23, # Assuming count was stable
            'count_jan_1_2023': 23,
            'has_disconnect_policy': True,
            'has_monitor_policy': True,
            'refuses_distribution': True,
        },
        {
            'name': 'C',
            'count_jan_1_2022': 1000,
            'count_jan_1_2023': 987, # Still > 25
            'has_disconnect_policy': True,
            'has_monitor_policy': True,
            'updates_monthly': True,
            'distributes_quarterly': True,
        },
        {
            'name': 'D',
            'count_jan_1_2022': 30,
            'count_jan_1_2023': 24, # 30 - 6
            'has_disconnect_policy': True,
            'has_monitor_policy': False,
        },
        {
            'name': 'E',
            'count_jan_1_2022': 22,
            'count_jan_1_2023': 18, # 22 - 2 - 2
            'has_disconnect_policy': False,
            'has_monitor_policy': False,
        }
    ]

    non_compliant_employers = []
    for emp in employers:
        is_compliant, name = analyze_compliance(emp)
        if not is_compliant:
            non_compliant_employers.append(name)
            
    print("--- FINAL RESULT ---")
    print(f"The analysis identified the following employer(s) as non-compliant: {non_compliant_employers}")
    print("\nEmployer D is clearly not in compliance because they failed to create a mandatory policy (Electronic Monitoring) by the deadline.")
    print("Employer C is also not in compliance due to a flawed policy update distribution schedule.")
    print("However, the failure to create a policy at all (Employer D) is the more fundamental and unambiguous violation.")
    final_answer = 'D'
    print(f"The most definitively non-compliant employer is {final_answer}.")
    
    # Required final output format
    print(f"\n<<<{final_answer}>>>")


# Execute the solution
solve_task()