import sys

# Define the legal requirements in Ontario as of the check date (Jan 2, 2023)
POLICY_EMPLOYEE_THRESHOLD = 25
HEADCOUNT_TRIGGER_DATE = "January 1, 2022"

def analyze_employer_compliance():
    """
    Analyzes each employer scenario based on Ontario's employment laws.
    """

    employers = [
        {
            "id": "A",
            "headcount_jan1_2022": 20,
            "has_disconnect_policy": False,
            "has_monitoring_policy": False,
            "distributes_policy": True,
            "notes": "Employed 20 people on Jan 1, 2022. Hired more in 2023."
        },
        {
            "id": "B",
            "headcount_jan1_2022": None,  # Not specified, current headcount is 23
            "has_disconnect_policy": True,
            "has_monitoring_policy": True,
            "distributes_policy": False,
            "notes": "Refuses to distribute policies to new employees. Headcount on the trigger date is ambiguous."
        },
        {
            "id": "C",
            "headcount_jan1_2022": 1000,
            "has_disconnect_policy": True,
            "has_monitoring_policy": True,
            "distributes_policy": True,
            "notes": "Employed 1,000 people and has implemented and distributed all required policies."
        },
        {
            "id": "D",
            "headcount_jan1_2022": 30,
            "has_disconnect_policy": True,
            "has_monitoring_policy": False,
            "distributes_policy": True,
            "notes": "Employed 30 people on Jan 1, 2022, but did not develop an electronic monitoring policy."
        },
        {
            "id": "E",
            "headcount_jan1_2022": 22,
            "has_disconnect_policy": False,
            "has_monitoring_policy": False,
            "distributes_policy": True,
            "notes": "Employed 22 people on Jan 1, 2022. Headcount changed later in the year."
        }
    ]

    non_compliant_employer = None

    print(f"Analyzing employer compliance based on the employee threshold of {POLICY_EMPLOYEE_THRESHOLD} on {HEADCOUNT_TRIGGER_DATE}.\n")

    for emp in employers:
        print(f"--- Analyzing Employer {emp['id']} ---")
        
        # Handle the ambiguous case B
        if emp['id'] == 'B':
            print("Note: The headcount on Jan 1, 2022, is not specified. If it were 25 or more, the refusal to distribute policies would be a violation. However, the requirement to have the policy is not definitively established by the facts provided.")
            print("Status: Potentially non-compliant, but ambiguous.\n")
            continue

        headcount = emp['headcount_jan1_2022']
        print(f"Headcount on {HEADCOUNT_TRIGGER_DATE}: {headcount}")
        
        is_required = headcount >= POLICY_EMPLOYEE_THRESHOLD
        print(f"Is a policy required? Checking if {headcount} >= {POLICY_EMPLOYEE_THRESHOLD}. Result: {is_required}")

        if not is_required:
            print("Conclusion: The employer is not required to have these policies as their headcount is below the threshold.")
            print("Status: Compliant\n")
            continue

        # If policies are required, check for their existence
        if is_required:
            compliant = True
            if not emp['has_disconnect_policy']:
                print("Violation: 'Disconnecting from Work' policy is required but not present.")
                compliant = False
            if not emp['has_monitoring_policy']:
                print("Violation: 'Electronic Monitoring' policy is required but not present.")
                compliant = False
            if not emp['distributes_policy']:
                print("Violation: Required policies are not being distributed.")
                compliant = False
            
            if compliant:
                print("Conclusion: All required policies are in place and distributed.")
                print("Status: Compliant\n")
            else:
                print("Conclusion: The employer failed to meet all legal requirements.")
                print("Status: NOT IN COMPLIANCE\n")
                if non_compliant_employer is None:
                     non_compliant_employer = emp['id']


    print("--- Final Result ---")
    if non_compliant_employer:
        print(f"Based on a clear and unambiguous violation, Employer {non_compliant_employer} is not in compliance with applicable employment laws.")
    else:
        print("All employers with clearly defined circumstances appear to be in compliance.")


if __name__ == '__main__':
    analyze_employer_compliance()