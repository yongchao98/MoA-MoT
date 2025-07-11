def analyze_bee_experiments():
    """
    Analyzes experimental data on honeybees and fungi to determine the correct conclusion.
    """

    # --- Data from the problem description ---
    # Mortality rates (%)
    mortality_non_infected = 10
    mortality_fungus_a = 35  # For most pollens
    mortality_fungus_b = 20
    mortality_fungus_c = 10

    # Productivity data for Fungus C (Eggs laid)
    productivity_c = {
        'buck': {'not_infected': 45, 'infected': 60},
        'lavender': {'not_infected': 30, 'infected': 50}
    }


    # --- Step-by-step Analysis ---

    print("Step 1: Determine if Fungus A is a pathogen.")
    print(f"The mortality rate of non-infected bees is {mortality_non_infected}%.")
    print(f"The mortality rate of bees infected with Fungus A is {mortality_fungus_a}%.")
    is_pathogen_a = mortality_fungus_a > mortality_non_infected
    print(f"Is Fungus A a pathogen? ({mortality_fungus_a} > {mortality_non_infected}) -> {is_pathogen_a}")
    if is_pathogen_a:
        print("Conclusion: Fungus A causes a significant increase in mortality, so it is a pathogen.\n")
    else:
        print("Conclusion: Fungus A does not increase mortality, so it is not a pathogen.\n")


    print("Step 2: Determine if Fungus B is a pathogen.")
    print(f"The mortality rate of non-infected bees is {mortality_non_infected}%.")
    print(f"The mortality rate of bees infected with Fungus B is {mortality_fungus_b}%.")
    is_pathogen_b = mortality_fungus_b > mortality_non_infected
    print(f"Is Fungus B a pathogen? ({mortality_fungus_b} > {mortality_non_infected}) -> {is_pathogen_b}")
    if is_pathogen_b:
        print("Conclusion: Fungus B causes an increase in mortality, so it is a pathogen.\n")
    else:
        print("Conclusion: Fungus B does not increase mortality, so it is not a pathogen.\n")


    print("Step 3: Determine the nature of Fungus C.")
    print(f"The mortality rate of non-infected bees is {mortality_non_infected}%.")
    print(f"The mortality rate of bees infected with Fungus C is {mortality_fungus_c}%.")
    is_pathogen_c = mortality_fungus_c > mortality_non_infected
    print(f"Does Fungus C increase mortality? ({mortality_fungus_c} > {mortality_non_infected}) -> {is_pathogen_c}")

    # Further analysis of Fungus C
    print("\nAnalyzing productivity for Fungus C:")
    eggs_not_infected_buck = productivity_c['buck']['not_infected']
    eggs_infected_buck = productivity_c['buck']['infected']
    print(f"With buck pollen, egg count changes from {eggs_not_infected_buck} (not infected) to {eggs_infected_buck} (infected).")

    print("\nConclusion for Fungus C: Since mortality does not increase ({0} == {1}) and productivity can even increase, Fungus C is not a pathogen. It is best described as a commensal (or potentially mutualistic) organism.".format(mortality_fungus_c, mortality_non_infected))


    # --- Final Conclusion ---
    print("\n" + "="*50)
    print("Final Summary:")
    print("Fungus A is a pathogen.")
    print("Fungus B is a pathogen.")
    print("Fungus C is a commensal.")
    print("This matches the statement in option I.")
    print("="*50)


analyze_bee_experiments()
print("<<<I>>>")