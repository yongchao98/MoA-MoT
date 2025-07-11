def analyze_compliance():
    """
    Analyzes employers based on Ontario's employment laws for workplace policies
    as of January 2, 2023.
    """
    EMPLOYEE_THRESHOLD = 25
    print(f"Ontario Employment Law Analysis as of January 2, 2023")
    print(f"The legal threshold for requiring 'Disconnect from Work' and 'Electronic Monitoring' policies is {EMPLOYEE_THRESHOLD} employees on January 1st of a given year.\n")

    # --- Case A ---
    emp_a_count = 20
    is_a_compliant = emp_a_count < EMPLOYEE_THRESHOLD
    print("--- Employer A Analysis ---")
    print(f"On Jan 1, 2022, Employer A had {emp_a_count} employees.")
    print(f"Equation: {emp_a_count} < {EMPLOYEE_THRESHOLD}")
    print(f"Result: Since the number of employees is less than {EMPLOYEE_THRESHOLD}, no policies were required. Employer A is IN COMPLIANCE.\n")

    # --- Case B ---
    emp_b_count = 23
    is_b_compliant = emp_b_count < EMPLOYEE_THRESHOLD
    print("--- Employer B Analysis ---")
    print(f"Employer B has {emp_b_count} employees.")
    print(f"Equation: {emp_b_count} < {EMPLOYEE_THRESHOLD}")
    print(f"Result: Since the number of employees is less than {EMPLOYEE_THRESHOLD}, the laws do not apply. Employer B is IN COMPLIANCE.\n")

    # --- Case C ---
    emp_c_count = 1000
    is_c_compliant = (emp_c_count >= EMPLOYEE_THRESHOLD) and True # Has both policies
    print("--- Employer C Analysis ---")
    print(f"On Jan 1, 2022, Employer C had {emp_c_count} employees.")
    print(f"Equation: {emp_c_count} >= {EMPLOYEE_THRESHOLD}")
    print(f"Result: The threshold was met. The employer created and distributed both required policies. Employer C is IN COMPLIANCE.\n")

    # --- Case D ---
    emp_d_count = 30
    has_disconnect_policy = True
    has_monitoring_policy = False
    is_d_compliant = (emp_d_count < EMPLOYEE_THRESHOLD) or (has_disconnect_policy and has_monitoring_policy)
    print("--- Employer D Analysis ---")
    print(f"On Jan 1, 2022, Employer D had {emp_d_count} employees.")
    print(f"Equation: {emp_d_count} >= {EMPLOYEE_THRESHOLD}")
    print(f"Result: The threshold was met. The employer was required to have BOTH policies.")
    print(f"Disconnecting from Work Policy: Implemented -> {has_disconnect_policy}")
    print(f"Electronic Monitoring Policy: Implemented -> {has_monitoring_policy}")
    print("Conclusion: The employer failed to implement the mandatory electronic monitoring policy. Employer D is NOT IN COMPLIANCE.\n")

    # --- Case E ---
    emp_e_initial = 22
    emp_e_resignations = 2
    emp_e_retirements = 2
    # The count on Jan 1, 2022 was 22, so no policies needed for 2022.
    # We check the count for Jan 1, 2023 to be thorough.
    emp_e_count_jan_1_2023 = emp_e_initial - emp_e_resignations - emp_e_retirements
    is_e_compliant = emp_e_initial < EMPLOYEE_THRESHOLD and emp_e_count_jan_1_2023 < EMPLOYEE_THRESHOLD
    print("--- Employer E Analysis ---")
    print(f"On Jan 1, 2022, Employer E had {emp_e_initial} employees, which is less than {EMPLOYEE_THRESHOLD}, so they were compliant for 2022.")
    print(f"To check for 2023, we calculate the employee count on Jan 1, 2023:")
    print(f"Equation: {emp_e_initial} - {emp_e_resignations} (resigned) - {emp_e_retirements} (retired) = {emp_e_count_jan_1_2023}")
    print(f"Result: On Jan 1, 2023, the employee count was {emp_e_count_jan_1_2023}, which is still less than {EMPLOYEE_THRESHOLD}. Employer E is IN COMPLIANCE.\n")
    
    # Final Answer
    if not is_d_compliant:
        final_answer = "D"
    else:
        final_answer = "Error in logic"

    print("-------------------------")
    print(f"The employer that is not in compliance is Employer {final_answer}.")
    print("-------------------------")

if __name__ == '__main__':
    analyze_compliance()