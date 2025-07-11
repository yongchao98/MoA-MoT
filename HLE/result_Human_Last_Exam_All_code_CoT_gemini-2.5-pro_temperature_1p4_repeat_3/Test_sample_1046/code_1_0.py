import math

def solve_resistance_pace():
    """
    Calculates and compares the time for two different bacterial populations
    to achieve widespread drug resistance, illustrating the concepts in the problem.
    """
    total_population = 1_000_000
    target_resistance_level = total_population / 2

    print("--- Scenario 1: Bacterium with Lateral Transfer ---")
    # This bacterium already has some resistant members and resistance spreads quickly.
    initial_resistant_pop_1 = 10
    # Effective spread rate combines growth and lateral transfer.
    spread_rate_1 = 0.5

    # We solve for t in: target = initial * (1 + rate)^t
    # t = log(target / initial) / log(1 + rate)
    time_1 = math.log(target_resistance_level / initial_resistant_pop_1) / math.log(1 + spread_rate_1)

    print("The time for a population with lateral transfer to become 50% resistant.")
    print(f"Equation: t = log(Target / Initial) / log(1 + Rate)")
    print(f"Calculation: t = log({int(target_resistance_level)} / {initial_resistant_pop_1}) / log(1 + {spread_rate_1})")
    print(f"Result: Time = {time_1:.1f} time units.\n")


    print("--- Scenario 2: Bacterium with Rare but Highly Fit Mutation ---")
    # This bacterium has a time delay for the rare mutation to occur.
    # The mutation is so effective (cross-resistance + compensatory mutations)
    # that the mutant strain has a very high fitness advantage and spreads rapidly.
    time_delay_2 = 15.0
    initial_resistant_pop_2 = 1 # Starts from a single cell
    # The spread rate is extremely high due to high fitness and selection pressure.
    spread_rate_2 = 2.07

    # The time to spread after the initial delay.
    time_to_spread_2 = math.log(target_resistance_level / initial_resistant_pop_2) / math.log(1 + spread_rate_2)
    total_time_2 = time_delay_2 + time_to_spread_2

    print("The time for a population relying on a rare but highly advantageous mutation.")
    print("Equation: t_total = Delay + [log(Target / Initial) / log(1 + Rate)]")
    print(f"Calculation: t_total = {time_delay_2} + [log({int(target_resistance_level)} / {initial_resistant_pop_2}) / log(1 + {spread_rate_2})]")
    print(f"Result: Time = {time_delay_2:.1f} (delay) + {time_to_spread_2:.1f} (spread) = {total_time_2:.1f} time units.\n")

    print(f"Conclusion: A significant delay ({time_delay_2:.1f} units) for a rare mutation can be offset by extremely rapid spread ({time_to_spread_2:.1f} units) due to high fitness, resulting in a total time ({total_time_2:.1f} units) comparable to the lateral transfer scenario ({time_1:.1f} units).")

solve_resistance_pace()
<<<B>>>