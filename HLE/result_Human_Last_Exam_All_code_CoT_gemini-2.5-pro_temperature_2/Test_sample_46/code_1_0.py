import collections

def analyze_rmbs_value(total_principal, recovery_rate, default_rate):
    """
    Calculates the impact of the default rate on the value of an RMBS pool.
    
    Args:
        total_principal (float): The total initial value of the mortgage pool.
        recovery_rate (float): The percentage of principal recovered from a defaulted loan after foreclosure.
        default_rate (float): The percentage of loans in the pool that default.
        
    Returns:
        A dictionary containing the calculated values.
    """
    
    defaulted_principal = total_principal * default_rate
    loss_per_dollar_of_default = 1 - recovery_rate
    total_loss = defaulted_principal * loss_per_dollar_of_default
    remaining_value = total_principal - total_loss
    
    # Using an ordered dict to maintain the print order
    results = collections.OrderedDict()
    results["Total Principal"] = total_principal
    results["Assumed Recovery Rate"] = recovery_rate
    results["Assumed Default Rate"] = default_rate
    results["Total Principal in Default"] = defaulted_principal
    results["Loss Given Default"] = total_loss
    results["Final RMBS Value"] = remaining_value
    
    return results

def print_scenario(name, results):
    """Helper function to print the results of a scenario."""
    print(f"--- {name} ---")
    for key, value in results.items():
        if key == "Final RMBS Value":
            print("-" * 20)
            print(f"{key:<25}: ${value:,.2f}")
        else:
            # Print rates as percentages
            if "Rate" in key:
                 print(f"{key:<25}: {value:.2%}")
            else:
                 print(f"{key:<25}: ${value:,.2f}")
    print("\n")


if __name__ == "__main__":
    # --- Parameters ---
    # Non-Agency RMBS are backed by mortgages. Their value depends on the homeowners paying their mortgages.
    # The value is primarily determined by: (E) Default rates and (G) Recovery rates.
    # Let's model this.
    PRINCIPAL = 100_000_000.00 # A hypothetical $100 million mortgage pool
    RECOVERY_RATE = 0.40      # A 40% recovery rate on foreclosed properties was common during the crisis.
    
    print("This simulation shows how the value of a non-agency RMBS is determined by its underlying default rate.")
    print("The key calculation is: Final Value = Total Principal - (Total Principal * Default Rate * (1 - Recovery Rate))\n")

    # --- Scenario 1: Pre-Crisis Expectations ---
    # Investors expected low default rates.
    normal_default_rate = 0.02 # 2% default rate
    normal_scenario_results = analyze_rmbs_value(PRINCIPAL, RECOVERY_RATE, normal_default_rate)
    print_scenario("Scenario 1: Low Default Rate (Expectation)", normal_scenario_results)
    
    # --- Scenario 2: The 2008 Crisis Reality ---
    # Subprime default rates surged to unprecedented levels.
    crisis_default_rate = 0.25 # 25% default rate
    crisis_scenario_results = analyze_rmbs_value(PRINCIPAL, RECOVERY_RATE, crisis_default_rate)
    print_scenario("Scenario 2: High Default Rate (Crisis Reality)", crisis_scenario_results)
    
    print("As you can see, the default rate is the most powerful variable in determining the security's final value.")
    print("While factors like FICO scores (C) and originator quality (F) are root causes of defaults,")
    print("the default rate itself is the direct measure of losses that 'determines the value'.")

<<<E>>>