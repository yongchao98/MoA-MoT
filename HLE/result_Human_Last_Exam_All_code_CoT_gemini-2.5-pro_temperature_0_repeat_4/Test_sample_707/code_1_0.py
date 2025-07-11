def calculate_flawed_probability():
    """
    This function replicates Agent A's flawed calculation based on a linear extrapolation.
    """
    # Initial parameters given in the problem
    p_initial = 0.01
    increase_per_quarter = 0.01
    years = 5
    quarters_per_year = 4

    # Agent A's calculation process
    total_quarters = years * quarters_per_year
    
    # The agent observes the increase after the first quarter.
    # He then extrapolates this for the remaining period.
    # There are 20 quarters in total. The calculation for the final probability
    # is based on the initial state plus 19 subsequent increases to reach p=0.2.
    num_increases = total_quarters - 1
    
    # The agent's flawed model: P_final = P_initial + num_increases * increase_per_quarter
    total_increase = num_increases * increase_per_quarter
    p_final = p_initial + total_increase

    print("Agent A's Flawed Calculation Demonstration:")
    print("The agent assumes a constant increase based on the first quarter's observation and extrapolates it linearly.")
    print("\n--- The Equation ---")
    print("Final Probability = Initial Probability + (Total Quarters - 1) * Per-Quarter Increase")
    
    # Outputting each number in the final equation as requested
    print("\n--- The Numbers ---")
    print(f"Final Probability = {p_initial} + ({total_quarters} - 1) * {increase_per_quarter}")
    print(f"Final Probability = {p_initial} + {num_increases} * {increase_per_quarter}")
    print(f"Final Probability = {p_initial} + {total_increase}")
    print(f"Final Probability = {p_final}")

calculate_flawed_probability()