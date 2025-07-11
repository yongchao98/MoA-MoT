def analyze_compliance():
    """
    Analyzes employer compliance with Ontario employment laws regarding specific policies.
    """
    employee_threshold = 25
    non_compliant_employer = None
    final_reasoning = ""

    # --- Employer A ---
    emp_a_count_jan1_2022 = 20
    is_a_compliant = emp_a_count_jan1_2022 < employee_threshold
    print("--- Analysis for Employer A ---")
    print(f"Employee count on Jan 1, 2022: {emp_a_count_jan1_2022}")
    print(f"Is count >= {employee_threshold}? {'Yes' if not is_a_compliant else 'No'}")
    print(f"Compliance Status: {'Compliant' if is_a_compliant else 'Not Compliant'}\n")

    # --- Employer B ---
    # The law applies if the employer has >= 25 employees. This employer has 23.
    emp_b_count = 23
    is_b_compliant = emp_b_count < employee_threshold
    print("--- Analysis for Employer B ---")
    print(f"Employee count: {emp_b_count}")
    print(f"Is count >= {employee_threshold}? {'Yes' if not is_b_compliant else 'No'}")
    print("Reason: The requirement to have and distribute policies is not triggered.")
    print(f"Compliance Status: {'Compliant' if is_b_compliant else 'Not Compliant'}\n")

    # --- Employer C ---
    emp_c_count_jan1_2022 = 1000
    # They are required to have policies and they do.
    is_c_compliant = True
    print("--- Analysis for Employer C ---")
    print(f"Employee count on Jan 1, 2022: {emp_c_count_jan1_2022}")
    print(f"Is count >= {employee_threshold}? {'Yes' if emp_c_count_jan1_2022 >= employee_threshold else 'No'}")
    print("Action: Has both required policies and distributes them.")
    print(f"Compliance Status: {'Compliant' if is_c_compliant else 'Not Compliant'}\n")

    # --- Employer D ---
    emp_d_count_jan1_2022 = 30
    has_monitoring_policy = False
    is_d_compliant = True
    print("--- Analysis for Employer D ---")
    print(f"Employee count on Jan 1, 2022: {emp_d_count_jan1_2022}")
    if emp_d_count_jan1_2022 >= employee_threshold:
        print(f"Is count >= {employee_threshold}? Yes")
        print("Requirement: Must have both 'disconnect from work' and 'electronic monitoring' policies.")
        if not has_monitoring_policy:
            is_d_compliant = False
            non_compliant_employer = "D"
            final_reasoning = (
                f"On January 1, 2022, the employer had {emp_d_count_jan1_2022} employees.\n"
                f"The legal threshold for requiring an electronic monitoring policy is {employee_threshold} employees.\n"
                f"Since {emp_d_count_jan1_2022} >= {employee_threshold}, the employer was required to have an electronic monitoring policy by October 11, 2022.\n"
                "The employer failed to create this policy and is therefore not in compliance."
            )
    else:
        print(f"Is count >= {employee_threshold}? No")
    print(f"Compliance Status: {'Compliant' if is_d_compliant else 'Not Compliant'}\n")

    # --- Employer E ---
    emp_e_count_jan1_2022 = 22
    is_e_compliant = emp_e_count_jan1_2022 < employee_threshold
    print("--- Analysis for Employer E ---")
    print(f"Employee count on Jan 1, 2022: {emp_e_count_jan1_2022}")
    print(f"Is count >= {employee_threshold}? {'Yes' if not is_e_compliant else 'No'}")
    print(f"Compliance Status: {'Compliant' if is_e_compliant else 'Not Compliant'}\n")

    # --- Final Answer ---
    print("="*30)
    print("FINAL CONCLUSION:")
    print(f"The employer not in compliance is: Employer {non_compliant_employer}")
    print("\nReasoning with numbers from the final equation:")
    print(final_reasoning)
    print("="*30)


if __name__ == "__main__":
    analyze_compliance()