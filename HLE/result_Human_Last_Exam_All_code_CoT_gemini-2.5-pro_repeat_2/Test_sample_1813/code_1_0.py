import math

def get_continued_fraction_coeffs(p, q):
    """
    Computes the continued fraction coefficients for the rational number p/q
    using the Euclidean algorithm.
    """
    if q == 0:
        raise ValueError("Denominator cannot be zero.")
    
    coeffs = []
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def compute_generalized_markov_cf():
    """
    Computes and prints the continued fraction associated with the 
    generalized Markov number m_{p/q}.
    """
    # Set the rational number for the problem
    p = 4
    q = 7

    # Step 1: Compute the continued fraction of p/q.
    # For 4/7, this results in [0, 1, 1, 3].
    coeffs = get_continued_fraction_coeffs(p, q)
    
    a0 = coeffs[0]
    W = coeffs[1:]

    # Step 2: Ensure the sequence W has an odd number of elements.
    # For W = (1, 1, 3), the length is 3 (odd), so no modification is needed.
    if len(W) % 2 == 0:
        # If the length were even, e.g., (a_1,...,a_n), we would transform it to
        # (a_1,...,a_n-1, 1) to make the length odd. This case is not needed for 4/7.
        an = W.pop()
        if an > 1:
            W.append(an - 1)
            W.append(1)
        else: # an == 1
             # This case requires a different rule, e.g., [..., a_{n-1}] -> [..., a_{n-1}+1]
             # but we'll stick to the simpler rule as it's not needed here.
             pass


    # Step 3: Construct the periodic part of the target continued fraction.
    # It is formed by concatenating W with its reverse, W^R.
    W_R = W[::-1]
    periodic_part = W + W_R

    # Step 4: Assemble and print the final continued fraction.
    # The final result is in the format [a_0; (periodic_part)], where parentheses
    # denote the repeating block. The "equation" for the number m_{4/7} is its
    # continued fraction representation, and we print all its constituent numbers.
    periodic_str = ", ".join(map(str, periodic_part))
    print(f"m_{p}/{q} = [{a0}; ({periodic_str})]")

compute_generalized_markov_cf()