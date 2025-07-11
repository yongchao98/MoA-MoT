import math

def get_continued_fraction_coefficients(p, q):
    """
    Computes the coefficients of the simple continued fraction for p/q.
    For a rational in (0,1), this is [a1, a2, ...].
    """
    # We expect a rational number p/q in (0, 1), so the first term is 0.
    # We are interested in the coefficients a_1, a_2, ...
    if q == 0:
        return []
    if p > q or p < 0:
        # This implementation is simplified for 0 < p/q < 1
        # The integer part is handled by the caller if necessary.
        pass

    coeffs = []
    temp_p, temp_q = p, q
    while temp_q != 0:
        a = temp_p // temp_q
        coeffs.append(a)
        temp_p, temp_q = temp_q, temp_p % temp_q
        
    # For p/q in (0,1), the first coefficient is 0. We return the rest.
    return coeffs[1:]

def compute_generalized_markov_number(coeffs):
    """
    Computes the generalized Markov number using the recurrence relation
    m_k = a_k * m_{k-1} - m_{k-2}.
    Returns the sequence of intermediate Markov numbers.
    """
    if not coeffs:
        # Base case for empty continued fraction
        return [1]

    # m_0 = 1 (by convention, as m of empty list)
    # m_1 = m([a_1]) = a_1
    markov_sequence = [1, coeffs[0]]

    for i in range(1, len(coeffs)):
        a_k = coeffs[i]
        m_k_minus_1 = markov_sequence[-1]
        m_k_minus_2 = markov_sequence[-2]
        
        m_k = a_k * m_k_minus_1 - m_k_minus_2
        markov_sequence.append(m_k)
        
    return markov_sequence

def solve():
    """
    Main function to compute the generalized Markov number for 4/7
    and print the details of the calculation.
    """
    p, q = 4, 7

    # Step 1: Find the continued fraction of 4/7.
    cf_coeffs = get_continued_fraction_coefficients(p, q)
    
    # Step 2: Compute the sequence of generalized Markov numbers.
    markov_sequence = compute_generalized_markov_number(cf_coeffs)

    # Step 3: Present the final calculation step as the "final equation".
    final_m = markov_sequence[-1]
    a_n = cf_coeffs[-1]
    m_n_minus_1 = markov_sequence[-2]
    m_n_minus_2 = markov_sequence[-3]

    print(f"The continued fraction for {p}/{q} is [0; {', '.join(map(str, cf_coeffs))}].")
    print(f"The sequence of generalized Markov numbers for the prefixes of the continued fraction is:")
    print(f"m(empty) = {markov_sequence[0]}")
    for i, c in enumerate(cf_coeffs):
        prefix = f"[{', '.join(map(str, cf_coeffs[:i+1]))}]"
        print(f"m({prefix}) = {markov_sequence[i+2]}")

    print("\nThe final equation computes the last term:")
    print(f"m([{', '.join(map(str, cf_coeffs))}]) = {a_n} * m([{', '.join(map(str, cf_coeffs[:-1]))}]) - m([{', '.join(map(str, cf_coeffs[:-2]))}])")
    print("Plugging in the numbers:")
    print(f"{final_m} = {a_n} * {m_n_minus_1} - {m_n_minus_2}")
    
    print("\nThe numbers in the final equation are:")
    print(final_m)
    print(a_n)
    print(m_n_minus_1)
    print(m_n_minus_2)


solve()