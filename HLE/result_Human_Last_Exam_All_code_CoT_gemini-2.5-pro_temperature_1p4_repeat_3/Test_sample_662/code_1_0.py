def calculate_vaccine_coverage_threshold():
    """
    Calculates the theoretical vaccine coverage threshold required to prevent an epidemic.
    """
    # Input parameters based on the user's information
    R0 = 3.0  # Basic reproduction number
    breakthrough_rate = 0.06  # 6% of vaccinated people get infected

    # Step 1: Calculate Vaccine Effectiveness (VE)
    # VE is the reduction in risk of infection for vaccinated individuals.
    vaccine_effectiveness = 1 - breakthrough_rate

    # Step 2: Calculate the critical vaccination coverage threshold (p_c)
    # The formula is p_c = (1 - 1/R0) / VE
    critical_coverage = (1 - 1 / R0) / vaccine_effectiveness

    # Step 3: Print the explanation and the final equation with numbers
    print("The formula to calculate the critical vaccination coverage (p_c) for a non-perfect vaccine is:")
    print("p_c = (1 - 1/R0) / VE\n")
    
    print("Where:")
    print(f"R0 (Basic Reproduction Number) = {R0}")
    print(f"VE (Vaccine Effectiveness) = 1 - Breakthrough Rate = {vaccine_effectiveness:.2f}\n")

    print("Plugging the numbers into the formula:")
    print(f"p_c = (1 - 1/{R0}) / {vaccine_effectiveness:.2f}")
    
    # Calculate and show the intermediate step for clarity
    herd_immunity_threshold = 1 - 1/R0
    print(f"p_c = {herd_immunity_threshold:.4f} / {vaccine_effectiveness:.2f}")
    
    print(f"p_c = {critical_coverage:.4f}\n")
    
    print(f"Conclusion: To prevent the spread of the virus, the theoretical threshold for vaccine coverage is {critical_coverage:.2%}.")

if __name__ == '__main__':
    calculate_vaccine_coverage_threshold()