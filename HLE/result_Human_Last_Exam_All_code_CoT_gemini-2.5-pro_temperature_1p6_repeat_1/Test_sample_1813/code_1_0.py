import math

def get_cf_coefficients(p, q):
    """
    Computes the standard continued fraction coefficients for a rational number p/q.
    The standard form is the one where the last coefficient is > 1.
    The Euclidean algorithm naturally produces this form.
    This function returns the coefficients after the integer part.
    """
    if q == 0:
        raise ValueError("Denominator cannot be zero.")
    if p < 0 or q < 0:
        # For simplicity, this implementation assumes positive p, q
        raise ValueError("p and q must be positive.")

    coeffs = []
    # We are interested in p/q in (0,1), so the integer part is 0.
    # The algorithm starts with q/p.
    num, den = q, p
    while den != 0:
        quotient = num // den
        remainder = num % den
        coeffs.append(quotient)
        num, den = den, remainder
    return coeffs

def compute_markov_cf(p, q):
    """
    Computes the continued fraction associated with the generalized Markov number m_{p/q}.
    """
    # Step 1: Find the standard continued fraction coefficients of p/q.
    # The coefficients C = (c_1, ..., c_m).
    c = get_cf_coefficients(p, q)
    print(f"The standard continued fraction of {p}/{q} is [0; {', '.join(map(str, c))}].")

    # Step 2: Form the palindromic sequence P = (C, C^R).
    c_rev = c[::-1]  # Reverse the list of coefficients
    periodic_part = c + c_rev
    
    return periodic_part

# The given rational number is 4/7
p_val = 4
q_val = 7

# Compute the periodic part of the associated continued fraction
periodic_coeffs = compute_markov_cf(p_val, q_val)

# Step 3: Output the final result.
# The final continued fraction is of the form [0; P] where P is the periodic part.
# The problem asks to output each number in the final equation.
# We interpret this as printing the numbers in the periodic part of the fraction.
print("\nThe continued fraction associated with the generalized Markov number m_{4/7} is:")
# We represent the repeating part with an overline in standard notation,
# e.g., [0; 1, 1, 3, 3, 1, 1]
# Here we just print the numbers of this repeating part.
print("Final Equation: [0; " + ", ".join(map(str, periodic_coeffs)) + "] (the part after the semicolon repeats)")
print("The numbers in the periodic part are:")
for number in periodic_coeffs:
    print(number)
