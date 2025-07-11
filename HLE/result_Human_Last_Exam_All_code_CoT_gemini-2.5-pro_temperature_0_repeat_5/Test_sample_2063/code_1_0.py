def solve_quantum_dna_probability():
    """
    Calculates the limiting probability P(n) as n -> infinity.

    The problem asks for the probability that for a random DNA sequence of length n
    with 8 bases, the condition (count of base x) mod 8 != x for all x in {0..7}.

    As n -> infinity, the probability that (count of base x) mod 8 = k approaches 1/8
    for any k. Thus, the probability of the failure condition for a single base x,
    P(count(x) mod 8 = x), is 1/8.

    The probability of success for a single base x is 1 - 1/8 = 7/8.

    Since there are 8 independent conditions to satisfy (one for each base x=0 to 7),
    the total probability of success is (7/8)^8.
    """

    # The base numbers for the equation
    base_success = 7
    base_total = 8
    power = 8

    # Calculate the numerator and denominator
    numerator = base_success ** power
    denominator = base_total ** power

    # Calculate the final probability
    probability = numerator / denominator

    # Print the results as requested
    print(f"The final equation for the limiting probability is: ({base_success}^{power}) / ({base_total}^{power})")
    print(f"Numerator: {base_success}^{power} = {numerator}")
    print(f"Denominator: {base_total}^{power} = {denominator}")
    print(f"The limiting value is: {numerator} / {denominator} = {probability}")

solve_quantum_dna_probability()