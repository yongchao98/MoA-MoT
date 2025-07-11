def explain_violation():
    """
    This function demonstrates how Critical-Level Views violate the
    Non-Elitism condition of Arrhenius's impossibility theorem.
    """

    # 1. Define the parameters for our model of a critical-level view.
    critical_level = 20  # The welfare level a new life must exceed to be a net positive.
    epsilon = 0.01         # A very small amount of welfare.
    large_N = 1000000      # Represents an arbitrarily large number of people.

    # 2. Define the Non-Elitism condition (in plain English).
    # Non-Elitism: There is NO "elite" welfare level such that one person at that level
    # is better than ANY number of people at a slightly lower level.

    # 3. We create two populations to test this condition.
    #    Population A has one person at a proposed "elite" level (just above critical).
    #    Population B has many people at a "lower" level (just below critical).
    welfare_A = critical_level + epsilon
    welfare_B = critical_level - epsilon

    # 4. Calculate the value of each population using the critical-level rule.
    # The rule: Value = sum of (individual welfare - critical_level)
    
    # Value for Population A: One person
    value_A = welfare_A - critical_level
    
    # Value for Population B: N people
    value_B = large_N * (welfare_B - critical_level)

    # 5. Print the results and the explanation.
    print("--- Testing Critical-Level Views against the Non-Elitism Condition ---")
    print(f"Let's assume a critical welfare level of: {critical_level}")
    print("\nWe will check if a welfare level just above this acts as an 'elite' level.")
    print("-" * 65)

    print("Population A: 1 person with welfare just ABOVE the critical level.")
    print(f"This person's welfare is: {welfare_A}")
    print("The value of this population is calculated as: (welfare) - (critical_level)")
    # Outputting the numbers in the final equation as requested
    print(f"Equation: ({welfare_A}) - ({critical_level}) = {value_A:.2f}")
    
    print("-" * 65)
    
    print(f"Population B: {large_N:,} people with welfare just BELOW the critical level.")
    print(f"Each person's welfare is: {welfare_B}")
    print("The value of this population is calculated as: N * ((welfare) - (critical_level))")
    # Outputting the numbers in the final equation as requested
    print(f"Equation: {large_N} * (({welfare_B}) - ({critical_level})) = {value_B:,.2f}")
    
    print("-" * 65)
    print("\nConclusion:")
    if value_A > value_B:
        print(f"Population A (Value: {value_A:.2f}) is considered better than Population B (Value: {value_B:,.2f}).")
        print("\nThis outcome holds true no matter how large Population B is. In fact, adding more people to Population B makes its value even more negative.")
        print("This means the welfare level just above the critical level acts as an 'elite' level.")
        print("A single life there is preferred over any number of lives slightly below it.")
        print("\nThis directly violates the 'Non-Elitism' condition.")

explain_violation()
<<<C>>>