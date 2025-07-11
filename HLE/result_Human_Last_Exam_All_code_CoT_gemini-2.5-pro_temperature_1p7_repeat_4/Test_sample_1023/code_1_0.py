def solve_extinction_rate_comparison():
    """
    This function explains the step-by-step derivation to find the ratio
    of morphospecies extinction rate to evolutionary species extinction rate.
    """
    
    # Define symbolic rates for explanation
    lambda_rate = "\u03BB"  # Speciation rate
    mu_rate = "\u03BC"      # True extinction rate
    alpha_rate = "\u03B1"    # Anagenesis rate

    print("Step 1: Define the extinction rate for each species concept.")
    print("Let's denote the fundamental rates as:")
    print(f"  {lambda_rate}: true speciation rate (cladogenesis)")
    print(f"  {mu_rate}: true biological extinction rate")
    print(f"  {alpha_rate}: anagenetic pseudoextinction rate")
    print("-" * 30)

    # Formula for evolutionary species extinction rate
    mu_e_formula = f"{mu_rate} + {lambda_rate}"
    print("For an Evolutionary Species, a taxon goes extinct upon true extinction or speciation.")
    print(f"Extinction Rate (mu_e) = {mu_e_formula}")
    print("-" * 30)

    # Formula for morphospecies extinction rate
    mu_m_formula = f"{mu_rate} + {alpha_rate} + 0.5 * {lambda_rate}"
    print("For a Morphospecies, a taxon goes extinct upon true extinction, anagenesis, or bifurcating speciation (50% of branching events).")
    print(f"Extinction Rate (mu_m) = {mu_m_formula}")
    print("-" * 30)

    print("Step 2: Apply the constraint that 'all processes occur at the same rates'.")
    print("We interpret this as the total origination rate of new taxa being equal for both concepts.")
    
    # Formula for evolutionary species origination rate
    O_e_formula = f"2 * {lambda_rate}"
    print(f"Origination Rate (O_e) = {O_e_formula} (a speciation event creates 2 new taxa)")

    # Formula for morphospecies origination rate
    O_m_formula = f"(1 * 0.5*{lambda_rate}) + (2 * 0.5*{lambda_rate}) + (1 * {alpha_rate}) = 1.5*{lambda_rate} + {alpha_rate}"
    print(f"Origination Rate (O_m) = {O_m_formula}")

    # Equating the rates to find the relationship
    print(f"\nSetting O_e = O_m gives: 2*{lambda_rate} = 1.5*{lambda_rate} + {alpha_rate}")
    alpha_relation = f"0.5 * {lambda_rate}"
    print(f"Solving for {alpha_rate}, we find: {alpha_rate} = {alpha_relation}")
    print("-" * 30)

    print("Step 3: Substitute this relationship back into the formula for mu_m.")
    mu_m_updated_formula = f"{mu_rate} + ({alpha_relation}) + 0.5 * {lambda_rate}"
    print(f"mu_m = {mu_m_updated_formula}")
    
    final_mu_m_formula = f"{mu_rate} + {lambda_rate}"
    print(f"Simplifying this gives: mu_m = {final_mu_m_formula}")
    print("-" * 30)

    print("Step 4: Calculate the final ratio of mu_m to mu_e.")
    print(f"We are comparing mu_m = {final_mu_m_formula} with mu_e = {mu_e_formula}.")
    print(f"The final equation for the ratio is: ({final_mu_m_formula}) / ({mu_e_formula})")

    # The numbers in the final equation are symbolic, so we calculate the numeric result
    final_factor = 1.0
    
    print(f"Since both expressions are identical, the ratio is 1.")
    print("\n-------------------------------------------")
    print(f"Final Answer: The multiplicative factor is {int(final_factor)}.")
    print("-------------------------------------------")

solve_extinction_rate_comparison()