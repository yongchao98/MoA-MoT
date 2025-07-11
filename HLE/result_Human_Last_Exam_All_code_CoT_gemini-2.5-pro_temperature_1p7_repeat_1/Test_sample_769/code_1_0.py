import sys

def solve():
    """
    Analyzes employer scenarios against Ontario's employment laws
    regarding 'disconnecting from work' and 'electronic monitoring' policies
    as of January 2, 2023.
    """

    EMPLOYEE_THRESHOLD = 25
    DETERMINATION_DATE = "January 1, 2022"

    employers = [
        {
            "id": "A",
            "employees_on_jan_1_2022": 20,
            "has_disconnect_policy": False,
            "has_monitoring_policy": False, # Not required to have one
            "description": "Had 20 employees on Jan 1, 2022. Hired more in 2023."
        },
        {
            "id": "B",
            # The prompt says the employer "has 23 employees". We interpret this as the count
            # on the relevant determination date, Jan 1, 2022, was also below the threshold.
            "employees_on_jan_1_2022": 23,
            "has_disconnect_policy": True,
            "has_monitoring_policy": True,
            "refuses_distribution": True,
            "description": "Has 23 employees and policies, but refuses to distribute them."
        },
        {
            "id": "C",
            "employees_on_jan_1_2022": 1000,
            "has_disconnect_policy": True,
            "has_monitoring_policy": True,
            "refuses_distribution": False,
            "description": "Has 1000 employees and all required policies, which are distributed."
        },
        {
            "id": "D",
            "employees_on_jan_1_2022": 30,
            "has_disconnect_policy": True,
            "has_monitoring_policy": False,
            "refuses_distribution": False,
            "description": "Had 30 employees on Jan 1, 2022, but lacks an electronic monitoring policy."
        },
        {
            "id": "E",
            "employees_on_jan_1_2022": 22,
            "has_disconnect_policy": False,
            "has_monitoring_policy": False, # Not required to have one
            "description": "Had 22 employees on Jan 1, 2022. Headcount changed later."
        }
    ]

    non_compliant_employer = None
    reasoning = ""

    print("Analyzing Employer Compliance as of January 2, 2023...\n")
    print(f"Rule: Policies are required if employee count on {DETERMINATION_DATE} was >= {EMPLOYEE_THRESHOLD}.\n")

    for emp in employers:
        is_compliant = True
        analysis = []
        
        num_employees = emp["employees_on_jan_1_2022"]
        analysis.append(f"Employer {emp['id']} had {num_employees} employees on {DETERMINATION_DATE}.")

        if num_employees < EMPLOYEE_THRESHOLD:
            analysis.append(f"Since {num_employees} < {EMPLOYEE_THRESHOLD}, policies were not required.")
            is_compliant = True
        else:
            analysis.append(f"Since {num_employees} >= {EMPLOYEE_THRESHOLD}, policies for disconnecting from work and electronic monitoring were required.")
            
            # Check for Disconnecting from Work Policy
            if not emp["has_disconnect_policy"]:
                analysis.append("-> FAIL: Did not have a 'disconnecting from work' policy.")
                is_compliant = False

            # Check for Electronic Monitoring Policy
            if not emp["has_monitoring_policy"]:
                analysis.append("-> FAIL: Did not have an 'electronic monitoring' policy.")
                is_compliant = False
            
            # Check for distribution (only for employer B's specific case)
            if emp.get("refuses_distribution"):
                analysis.append("-> NOTE: Refusing to distribute a required policy is a violation.")
                # This logic is secondary, as the headcount makes it non-required.
                # If they were required, this would be a clear fail.
                
        status = "Compliant"
        if not is_compliant:
            status = "NOT COMPLIANT"
            non_compliant_employer = emp['id']
            reasoning = "\n".join(analysis)

        print(f"--- Analysis for Employer {emp['id']} ---")
        for line in analysis:
            print(line)
        print(f"Result: {status}\n")

    if non_compliant_employer:
        print("--- Conclusion ---")
        print(f"The employer not in compliance is {non_compliant_employer}.")
        print("Final Reasoning:")
        print(reasoning)

        # Output the answer in the specified format
        sys.stdout.write(f'<<<{non_compliant_employer}>>>')
    else:
        print("All employers appear to be compliant based on the analysis.")

solve()