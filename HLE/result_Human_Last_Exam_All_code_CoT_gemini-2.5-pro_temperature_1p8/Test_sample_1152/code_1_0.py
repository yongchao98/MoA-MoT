def analyze_bee_experiments():
    """
    Analyzes experimental data on honeybees, fungi, and pollen to determine the correct statement.
    """
    print("Step 1: Analyzing the nature of Fungus A.")
    mortality_rate_fungus_a_infected = 35  # percent
    mortality_rate_non_infected = 10      # percent
    
    print(f"The mortality rate for non-infected bees is {mortality_rate_non_infected}%.")
    print(f"With Fungus A infection (on most pollens), the mortality rate is {mortality_rate_fungus_a_infected}%.")
    
    is_pathogen_A = mortality_rate_fungus_a_infected > mortality_rate_non_infected
    if is_pathogen_A:
        print(f"Since {mortality_rate_fungus_a_infected} > {mortality_rate_non_infected}, Fungus A increases mortality and is therefore a pathogen.\n")
    else:
        print(f"Since {mortality_rate_fungus_a_infected} <= {mortality_rate_non_infected}, Fungus A is not a pathogen.\n")

    print("Step 2: Analyzing the nature of Fungus B.")
    mortality_rate_fungus_b_infected = 20  # percent
    
    print(f"The mortality rate for non-infected bees is {mortality_rate_non_infected}%.")
    print(f"With Fungus B infection, the mortality rate is {mortality_rate_fungus_b_infected}%.")
    
    is_pathogen_B = mortality_rate_fungus_b_infected > mortality_rate_non_infected
    if is_pathogen_B:
        print(f"Since {mortality_rate_fungus_b_infected} > {mortality_rate_non_infected}, Fungus B increases mortality and is therefore a pathogen.\n")
    else:
        print(f"Since {mortality_rate_fungus_b_infected} <= {mortality_rate_non_infected}, Fungus B is not a pathogen.\n")

    print("Step 3: Analyzing the nature of Fungus C.")
    mortality_rate_fungus_c_infected = 10  # percent
    
    print(f"The mortality rate for non-infected bees is {mortality_rate_non_infected}%.")
    print(f"With Fungus C infection, the mortality rate is {mortality_rate_fungus_c_infected}%.")

    is_pathogen_C = mortality_rate_fungus_c_infected > mortality_rate_non_infected
    if not is_pathogen_C:
        # A commensal is an organism that lives with another without harming it.
        print(f"Since the mortality rate {mortality_rate_fungus_c_infected}% is the same as the non-infected rate {mortality_rate_non_infected}%, Fungus C does not cause harm in terms of mortality.")
        print("This means Fungus C is not a pathogen. It can be classified as a commensal (or even symbiotic, as it increases productivity in some cases).\n")
    else:
        print("Fungus C is a pathogen.\n")

    print("Step 4: Final Conclusion.")
    print("Based on the analysis:")
    print("- Fungus A is a pathogen.")
    print("- Fungus B is a pathogen.")
    print("- Fungus C is a commensal.")
    print("This directly corresponds to answer choice I.")

analyze_bee_experiments()