def analyze_compliance_scenarios():
    """
    Analyzes several employer scenarios against Ontario employment law as of Jan 2, 2023.
    The code explicitly shows the numbers and comparisons used in the legal analysis.
    """
    employee_threshold = 25
    print(f"--- Legal Framework ---")
    print(f"Assessment Date: January 2, 2023")
    print(f"Rule: Employers with {employee_threshold} or more employees on January 1 of a year must have written policies on 'disconnecting from work' and 'electronic monitoring'.")
    print(f"For 2022, the deadlines were June 2 (disconnecting) and October 11 (monitoring).\n")

    # --- Scenario A ---
    print("--- Analysis of Employer A ---")
    employees_a = 20
    print(f"Employer A had {employees_a} employees on January 1, 2022.")
    is_required_a = employees_a >= employee_threshold
    print(f"Is a policy required? ({employees_a} >= {employee_threshold}) -> {is_required_a}")
    print("Conclusion: With fewer than 25 employees, Employer A was not required to have these policies. They are in compliance.\n")

    # --- Scenario B ---
    print("--- Analysis of Employer B ---")
    employees_b = 23
    print(f"Employer B has {employees_b} employees.")
    is_required_b = employees_b >= employee_threshold
    print(f"Is a policy required? ({employees_b} >= {employee_threshold}) -> {is_required_b}")
    print("Conclusion: With fewer than 25 employees, Employer B is not legally required to have or distribute these policies. They are in compliance.\n")

    # --- Scenario C ---
    print("--- Analysis of Employer C ---")
    employees_c = 1000
    print(f"Employer C had {employees_c} employees on January 1, 2022.")
    is_required_c = employees_c >= employee_threshold
    print(f"Is a policy required? ({employees_c} >= {employee_threshold}) -> {is_required_c}")
    print("Obligation: Required to have both policies and distribute changes within 30 days.")
    print("Violation: The employer updates policies monthly but distributes them quarterly (approx. every 90 days). A 90-day distribution schedule for monthly updates violates the 30-day rule.")
    print(f"Is the distribution timely? (90 days <= 30 days) -> {90 <= 30}")
    print("Conclusion: Employer C is NOT in compliance due to its policy update distribution schedule.\n")

    # --- Scenario D ---
    print("--- Analysis of Employer D ---")
    employees_d = 30
    fired_employees = 6
    remaining_employees = employees_d - fired_employees
    print(f"Employer D had {employees_d} employees on January 1, 2022.")
    is_required_d = employees_d >= employee_threshold
    print(f"Is a policy required? ({employees_d} >= {employee_threshold}) -> {is_required_d}")
    print(f"Note: Firing {fired_employees} employees later in the year does not remove the obligation for 2022. The count on Jan 1 is what matters.")
    print("Obligation: Required to have both policies by the deadlines in 2022.")
    print("Violation: The employer 'has not developed a policy on electronic monitoring'.")
    print("Conclusion: Employer D is NOT in compliance because it failed to create a legally required policy.\n")

    # --- Scenario E ---
    print("--- Analysis of Employer E ---")
    employees_e = 22
    print(f"Employer E had {employees_e} employees on January 1, 2022.")
    is_required_e = employees_e >= employee_threshold
    print(f"Is a policy required for 2022? ({employees_e} >= {employee_threshold}) -> {is_required_e}")
    employees_on_jan_1_2023 = 22 - 2 - 2
    is_required_2023 = employees_on_jan_1_2023 >= employee_threshold
    print(f"Employee count on Jan 1, 2023 was {employees_on_jan_1_2023}, so no policy is required for 2023 either.")
    print("Conclusion: Employer E has consistently been below the 25-employee threshold and is in compliance.\n")
    
    print("--- Final Determination ---")
    print("Both employers C and D are not in compliance.")
    print("However, Employer D's violation is a complete failure to create a required policy.")
    print("Employer C's violation is procedural (improper distribution of updates).")
    print("The most direct and significant act of non-compliance is that of Employer D.")

analyze_compliance_scenarios()