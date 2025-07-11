def calculate_involutions():
    """
    Calculates the number of involutions for several finite groups
    and compares them according to the given choices.
    """

    # --- Calculations for each group ---

    # PSL(3,4): n=3, q=4 (char 2)
    # Formula: q * (q + 1) * (q^3 - 1)
    q_psl34 = 4
    n_psl34 = q_psl34 * (q_psl34 + 1) * (q_psl34**3 - 1)

    # PSU(3,3): n=3, q=3
    # Formula: q^2 * (q^2 - q + 1)
    q_psu33 = 3
    n_psu33 = q_psu33**2 * (q_psu33**2 - q_psu33 + 1)

    # PSL(3,9): n=3, q=9 (char 3)
    # Formula: q^2 * (q^3 - 1)
    q_psl39 = 9
    n_psl39 = q_psl39**2 * (q_psl39**3 - 1)

    # PSL(4,3): Two classes of involutions with sizes 5280 and 2112.
    n_psl43 = 5280 + 2112

    # PSU(4,4): n=4, q=4 (char 2), so m=2
    # Formula: q^(2m-1)*(q^(2m)-1)/(q+1) + q^(2m-2)*(q^m+(-1)^m)*(q^(m-1)-(-1)^(m-1))
    q_psu44 = 4
    m_psu44 = 2
    term1 = q_psu44**(2 * m_psu44 - 1) * (q_psu44**(2 * m_psu44) - 1) // (q_psu44 + 1)
    term2 = q_psu44**(2 * m_psu44 - 2) * (q_psu44**m_psu44 + (-1)**m_psu44) * (q_psu44**(m_psu44 - 1) - (-1)**(m_psu44 - 1))
    n_psu44 = term1 + term2

    results = {
        "PSL(3,4)": n_psl34,
        "PSU(3,3)": n_psu33,
        "PSL(3,9)": n_psl39,
        "PSL(4,3)": n_psl43,
        "PSU(4,4)": n_psu44
    }

    print("Number of involutions in each group:")
    for group, num in results.items():
        print(f"- {group}: {num}")
    print("\nComparing the pairs from the answer choices:")

    # --- Comparisons for answer choices ---

    # Choice A
    group1_A, group2_A = "PSL(3,4)", "PSU(3,3)"
    val1_A, val2_A = results[group1_A], results[group2_A]
    equal_A = "==" if val1_A == val2_A else "!="
    print(f"A. {group1_A} vs {group2_A}: {val1_A} {equal_A} {val2_A}")

    # Choice B
    group1_B, group2_B = "PSL(3,9)", "PSL(4,3)"
    val1_B, val2_B = results[group1_B], results[group2_B]
    equal_B = "==" if val1_B == val2_B else "!="
    print(f"B. {group1_B} vs {group2_B}: {val1_B} {equal_B} {val2_B}")

    # Choice C
    group1_C, group2_C = "PSL(3,9)", "PSU(4,4)"
    val1_C, val2_C = results[group1_C], results[group2_C]
    equal_C = "==" if val1_C == val2_C else "!="
    print(f"C. {group1_C} vs {group2_C}: {val1_C} {equal_C} {val2_C}")

    # Choice D
    group1_D, group2_D = "PSL(3,4)", "PSL(3,9)"
    val1_D, val2_D = results[group1_D], results[group2_D]
    equal_D = "==" if val1_D == val2_D else "!="
    print(f"D. {group1_D} vs {group2_D}: {val1_D} {equal_D} {val2_D}")

    if val1_A == val2_A or val1_B == val2_B or val1_C == val2_C or val1_D == val2_D:
      pass # A correct choice would be handled here
    else:
      print("\nConclusion: None of the pairs have an equal number of involutions.")

calculate_involutions()