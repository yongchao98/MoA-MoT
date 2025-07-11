import sympy

def solve_particle_annihilation():
    """
    This function calculates the expected time until the second collision
    for a system of five annihilating random walks, as described in the problem.
    """
    # Define symbolic variables for the initial separations
    N1, M1, N2, M2 = sympy.symbols('N1 M1 N2 M2')

    # These represent the initial gaps between adjacent particles.
    # d1=p2-p1, d2=p3-p2, d3=p4-p3, d4=p5-p4
    gaps = [N1, M1, N2, M2]
    n = 5 # Initial number of particles

    # --- Step 1: Calculate the expected rate-1 time for the FIRST collision ---
    # Formula: E[S_1] = (1/(n-2)) * sum_{k=1}^{n-2} d_k * d_{k+1}
    # For n=5, this is E[S_1] = (1/3) * (d1*d2 + d2*d3 + d3*d4).
    e_s1_sum_term = 0
    for k in range(n - 2): # k from 0 to 2
        e_s1_sum_term += gaps[k] * gaps[k+1]
    E_S1 = e_s1_sum_term / (n - 2)

    # --- Step 2: Calculate the expected rate-1 time to go from 5 particles to 1 ---
    # This is the sum of the rate-1 times for the first and second collisions, E[S_total] = E[S1 + S2].
    # Formula: E[S_total] = sum_{1 <= i < j <= n-1} d_i * d_j
    E_S_total = 0
    num_gaps = n - 1
    for i in range(num_gaps):
        for j in range(i + 1, num_gaps):
            E_S_total += gaps[i] * gaps[j]

    # --- Step 3: Calculate the expected physical time E[tau] ---
    # The physical time tau = T1 + T2.
    # T1 (rate 1) = S1. T2 (rate 2) = S2/2. So tau = S1 + S2/2.
    # E[tau] = E[S1] + E[S2]/2
    # We also know E[S_total] = E[S1] + E[S2], so E[S2] = E[S_total] - E[S1].
    # Substituting: E[tau] = E[S1] + (E[S_total] - E[S1])/2 = (E[S1] + E[S_total])/2.
    E_tau = (E_S1 + E_S_total) / 2

    # --- Step 4: Expand the expression for clarity and print the result ---
    final_expression = sympy.expand(E_tau)

    print("The expected time for the second collision is E[tau] =", end="")
    # The .as_ordered_terms() method provides the terms in a canonical order for consistent output.
    terms = final_expression.as_ordered_terms()
    for i, term in enumerate(terms):
        if i > 0 and term.could_extract_minus_sign():
            print(" -", end="")
            term = -term
        else:
            print(" +", end="")
        
        # A term can be split into its coefficient and variable parts
        coeff, factors = term.as_coeff_Mul()
        
        # Format the coefficient as a fraction, e.g., "1/2", "2/3"
        if isinstance(coeff, sympy.core.numbers.Half):
            coeff_str = "1/2"
        elif isinstance(coeff, sympy.core.numbers.Rational):
            coeff_str = f"{coeff.p}/{coeff.q}"
        else:
            coeff_str = str(coeff)

        # Format the variable part by joining factors with '*'
        vars_str = "*".join(map(str, factors.as_ordered_factors()))

        print(f" ({coeff_str})*{vars_str}", end="")
    print() # Final newline

if __name__ == '__main__':
    solve_particle_annihilation()