import datetime

def analyze_employers():
    """
    Analyzes several employer scenarios to determine compliance with
    Ontario employment laws as of January 2, 2023.
    """
    assessment_date = datetime.date(2023, 1, 2)
    employee_threshold = 25

    print(f"--- Analysis based on Ontario employment laws as of {assessment_date} ---\n")
    print(f"Key Rule: Policies on 'Disconnecting from Work' and 'Electronic Monitoring' are required if an employer has {employee_threshold} or more employees on January 1st of the year.\n")

    non_compliant_employer = None

    # Scenario A
    print("--- Analyzing Employer A ---")
    employees_jan_1_2022_a = 20
    print(f"Employer A had {employees_jan_1_2022_a} employees on January 1, 2022.")
    if employees_jan_1_2022_a < employee_threshold:
        print(f"Since {employees_jan_1_2022_a} is less than {employee_threshold}, Employer A was NOT required to have a 'Disconnecting from Work' or 'Electronic Monitoring' policy in 2022.")
        # Check Jan 1, 2023 count
        employees_jan_1_2023_a = 20
        print(f"The employee count on January 1, 2023 was {employees_jan_1_2023_a}, which is also below the threshold.")
        print("Result for A: IN COMPLIANCE.\n")
    else:
        print("Result for A: ANALYSIS ERROR.\n")

    # Scenario B
    print("--- Analyzing Employer B ---")
    current_employees_b = 23
    print(f"Employer B currently has {current_employees_b} employees.")
    print("The employee count on January 1, 2022, is not specified. However, their current count is below the threshold.")
    if current_employees_b < employee_threshold:
        print(f"If their count was also below {employee_threshold} on Jan 1, 2022, the laws requiring the policies do not apply to them.")
        print("Therefore, failure to distribute a policy they are not legally required to have does not constitute non-compliance with these specific laws.")
        print("Result for B: ASSUMED IN COMPLIANCE.\n")
    else:
        print("Result for B: NON-COMPLIANT (based on assumed prior headcount), but data is incomplete.\n")

    # Scenario C
    print("--- Analyzing Employer C ---")
    employees_jan_1_2022_c = 1000
    print(f"Employer C had {employees_jan_1_2022_c} employees on January 1, 2022.")
    if employees_jan_1_2022_c >= employee_threshold:
        print(f"Since {employees_jan_1_2022_c} is greater than or equal to {employee_threshold}, both policies are required.")
        print("The employer has written and appropriately distributed both policies. They go above and beyond legal requirements.")
        print("Result for C: IN COMPLIANCE.\n")
    else:
        print("Result for C: ANALYSIS ERROR.\n")


    # Scenario D
    print("--- Analyzing Employer D ---")
    employees_jan_1_2022_d = 30
    fired_employees_d = 6
    final_employees_2022_d = employees_jan_1_2022_d - fired_employees_d
    print(f"Employer D had {employees_jan_1_2022_d} employees on January 1, 2022.")
    print(f"The trigger for the policy requirement is based on the Jan 1 count ({employees_jan_1_2022_d}), not the later count of {final_employees_2022_d}.")

    if employees_jan_1_2022_d >= employee_threshold:
        print(f"Since {employees_jan_1_2022_d} is greater than or equal to {employee_threshold}, both policies were required.")
        print("The employer has a 'Disconnecting from Work' policy. This is compliant.")
        print("However, the employer 'has not developed a policy on electronic monitoring'.")
        print("The law required the employer to HAVE this policy in place by October 11, 2022.")
        print(f"By the assessment date of {assessment_date}, this policy is missing.")
        print("Result for D: NOT IN COMPLIANCE.\n")
        non_compliant_employer = 'D'
    else:
        print("Result for D: ANALYSIS ERROR.\n")

    # Scenario E
    print("--- Analyzing Employer E ---")
    employees_jan_1_2022_e = 22
    departed_2022_e = 4
    employees_jan_1_2023_e = employees_jan_1_2022_e - departed_2022_e
    print(f"Employer E had {employees_jan_1_2022_e} employees on January 1, 2022.")
    if employees_jan_1_2022_e < employee_threshold:
        print(f"Since {employees_jan_1_2022_e} is less than {employee_threshold}, Employer E was NOT required to have a policy in 2022.")
        print(f"Their employee count on Jan 1, 2023 was {employees_jan_1_2023_e}, which is also below the threshold.")
        print("Result for E: IN COMPLIANCE.\n")
    else:
        print("Result for E: ANALYSIS ERROR.\n")
        
    print(f"--- Conclusion ---")
    if non_compliant_employer:
        print(f"The analysis shows that Employer {non_compliant_employer} is not in compliance.")
    else:
        print("Could not definitively determine a single non-compliant employer.")


if __name__ == '__main__':
    analyze_employers()