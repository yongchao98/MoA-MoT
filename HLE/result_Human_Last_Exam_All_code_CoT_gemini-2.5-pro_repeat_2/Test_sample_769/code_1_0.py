import sys

def analyze_employer_compliance():
    """
    Analyzes employer compliance based on Ontario's laws for work-related policies.
    """
    # --- Legal Parameters ---
    EMPLOYEE_THRESHOLD = 25
    CRITICAL_DATE_2022 = "January 1, 2022"
    EVALUATION_DATE = "January 2, 2023"

    # --- Employer Data ---
    # We deduce the state of each employer based on the provided text.
    employers = [
        {
            "id": "A",
            "description": "Employed 20 people on Jan 1, 2022. Hired more in 2023.",
            "employees_on_crit_date": 20,
            "has_disconnect_policy": False,
            "has_monitoring_policy": False,
        },
        {
            "id": "B",
            "description": "Has 23 employees. Has policies but refuses to distribute to new hires.",
            "employees_on_crit_date": 23,
            "has_disconnect_policy": True,
            "has_monitoring_policy": True,
            "distributed_policy": False,
        },
        {
            "id": "C",
            "description": "Has over 987 employees. Has and distributes both policies.",
            "employees_on_crit_date": 1000,
            "has_disconnect_policy": True,
            "has_monitoring_policy": True,
            "distributed_policy": True,
        },
        {
            "id": "D",
            "description": "Employed 30 people on Jan 1, 2022, but not all policies are in place.",
            "employees_on_crit_date": 30,
            "has_disconnect_policy": True,
            "has_monitoring_policy": False,
        },
        {
            "id": "E",
            "description": "Had 22 employees on Jan 1, 2022. Headcount later decreased.",
            "employees_on_crit_date": 22,
            "has_disconnect_policy": False,
            "has_monitoring_policy": False,
        },
    ]

    non_compliant_employer = None

    print(f"Analyzing compliance as of {EVALUATION_DATE} based on the employee count on {CRITICAL_DATE_2022}.\n")
    print(f"The legal threshold for requiring policies is {EMPLOYEE_THRESHOLD} employees.\n")
    
    for emp in employers:
        print(f"--- Checking Employer {emp['id']} ---")
        print(f"Description: {emp['description']}")
        
        count = emp['employees_on_crit_date']
        is_required = count >= EMPLOYEE_THRESHOLD

        print(f"Employee count on {CRITICAL_DATE_2022}: {count}")
        print(f"Is a policy required? ({count} >= {EMPLOYEE_THRESHOLD}): {is_required}")

        if not is_required:
            print("Result: COMPLIANT. The employer was below the threshold, so no policies were legally required.\n")
            continue

        # If policies were required, check if they exist.
        disconnect_compliant = emp['has_disconnect_policy']
        monitoring_compliant = emp['has_monitoring_policy']
        
        print(f"Has 'Disconnecting from Work' policy: {disconnect_compliant}")
        print(f"Has 'Electronic Monitoring' policy: {monitoring_compliant}")
        
        if disconnect_compliant and monitoring_compliant:
            print("Result: COMPLIANT. The employer met the threshold and has both required policies.\n")
        else:
            print("Result: NOT COMPLIANT. The employer met the threshold but is missing one or more required policies.\n")
            non_compliant_employer = emp

    if non_compliant_employer:
        print("--- FINAL CONCLUSION ---")
        emp_id = non_compliant_employer['id']
        emp_count = non_compliant_employer['employees_on_crit_date']
        
        print(f"Employer {emp_id} is not in compliance.")
        print(f"The reason is that their employee count on {CRITICAL_DATE_2022} was {emp_count}.")
        # The final equation as requested
        print(f"The compliance check is: {emp_count} >= {EMPLOYEE_THRESHOLD}, which is True.")
        print(f"Because the threshold was met, the employer was required to have an electronic monitoring policy by Oct 11, 2022, but they failed to create one.")

        # This will be captured by the testing framework
        sys.stdout.flush()
        print(f"<<<{emp_id}>>>")

if __name__ == "__main__":
    analyze_employer_compliance()