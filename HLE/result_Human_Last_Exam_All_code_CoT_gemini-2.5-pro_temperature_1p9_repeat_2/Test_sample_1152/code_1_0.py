# Data from the experiments
control_mortality_rate = 10  # in percent

# Fungus A data (e.g., with buck pollen)
fungus_a_mortality_rate = 35 # in percent

# Fungus B data (e.g., with buck pollen)
fungus_b_mortality_rate = 20 # in percent

# Fungus C data (e.g., with buck pollen)
fungus_c_mortality_rate = 10 # in percent
fungus_c_productivity_infected = 60 # eggs
fungus_c_productivity_uninfected = 45 # eggs


def analyze_statement_I():
    """
    This function analyzes the claims made in statement I.
    I: Fungus A and B are pathogens. Fungus C is a commensal.
    """
    print("Analyzing Statement I: 'Fungus A and B are pathogens. Fungus C is a commensal.'\n")

    # Part 1: Fungus A is a pathogen
    print("Part 1: Is Fungus A a pathogen?")
    is_pathogen_A = fungus_a_mortality_rate > control_mortality_rate
    print(f"A pathogen increases the mortality rate above the control.")
    print(f"The mortality rate with Fungus A is {fungus_a_mortality_rate}%, which is greater than the control mortality rate of {control_mortality_rate}%.")
    print(f"Equation: {fungus_a_mortality_rate} > {control_mortality_rate} is {is_pathogen_A}.")
    print("Conclusion: Fungus A is a pathogen.\n")


    # Part 2: Fungus B is a pathogen
    print("Part 2: Is Fungus B a pathogen?")
    is_pathogen_B = fungus_b_mortality_rate > control_mortality_rate
    print(f"The mortality rate with Fungus B is {fungus_b_mortality_rate}%, which is greater than the control mortality rate of {control_mortality_rate}%.")
    print(f"Equation: {fungus_b_mortality_rate} > {control_mortality_rate} is {is_pathogen_B}.")
    print("Conclusion: Fungus B is a pathogen.\n")


    # Part 3: Fungus C is a commensal
    print("Part 3: Is Fungus C a commensal?")
    is_not_harmful_C = fungus_c_mortality_rate <= control_mortality_rate
    print(f"A commensal does not harm the host, meaning the mortality rate is not increased.")
    print(f"The mortality rate with Fungus C is {fungus_c_mortality_rate}%, which is equal to the control mortality rate of {control_mortality_rate}%.")
    print(f"Equation: {fungus_c_mortality_rate} <= {control_mortality_rate} is {is_not_harmful_C}.")
    print(f"Additionally, productivity in some cases increases (e.g., from {fungus_c_productivity_uninfected} to {fungus_c_productivity_infected} eggs), suggesting no harm and possible benefit.")
    print("Conclusion: Fungus C is a commensal.\n")

    # Final Conclusion
    if is_pathogen_A and is_pathogen_B and is_not_harmful_C:
        print("Final Verdict: Statement I is fully supported by the experimental data.")

analyze_statement_I()
<<<I>>>