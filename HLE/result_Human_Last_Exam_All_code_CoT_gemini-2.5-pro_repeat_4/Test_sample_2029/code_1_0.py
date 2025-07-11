import math

def safe_comb(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k > n or k < 0, which is standard.
    Handles the case where n < 0, which math.comb does not.
    In our formula, the 'n' part '61 - 7*i - 9*j' can be less than 'k' (13).
    math.comb(n,k) returns 0 if k>n, so we only need to guard against n<0.
    """
    if k < 0 or n < k:
        return 0
    return math.comb(n, k)

def solve_sequences():
    """
    Calculates the total number of different sequences (a_1,...,a_8, b_1,...,b_6)
    that can occur.
    """
    total_sequences = 0
    m = 8  # members in team A
    n = 6  # members in team B
    
    # The formula is derived from finding the coefficient of x^48 in the expansion of
    # (1+...+x^6)^8 * (1+...+x^8)^6.
    # Formula: sum_{i=0 to 8} sum_{j=0 to 6} C(8,i) * C(6,j) * (-1)**(i+j) * C(61 - 7*i - 9*j, 13)
    
    # The equation for the final answer is built by summing up all terms.
    # We will compute this sum.
    
    # Loop over the first binomial expansion term for (1-x^7)^8
    for i in range(m + 1):
        # Loop over the second binomial expansion term for (1-x^9)^6
        for j in range(n + 1):
            
            c_m_i = math.comb(m, i)
            c_n_j = math.comb(n, j)
            
            sign = (-1)**(i + j)
            
            # This is the 'n' for C(n,k) where k=13
            # from the expansion of (1-x)^(-14)
            # The power of x we need is 48 - 7*i - 9*j
            # The coefficient of x^k in (1-x)^-d is C(k+d-1, d-1)
            # k = 48 - 7i - 9j, d = 14
            # n_val = k + d - 1 = 48 - 7i - 9j + 14 - 1 = 61 - 7i - 9j
            n_val = m * n - (n + 1) * i - (m + 1) * j + (m + n - 1)
            k_val = m + n - 1
            
            c_term = safe_comb(n_val, k_val)
            
            term = c_m_i * c_n_j * sign * c_term
            total_sequences += term
            
    # The final equation is the result of the summation.
    # We print it in a descriptive way.
    final_equation = f"Total number of possible sequences = {total_sequences}"
    print(final_equation)

solve_sequences()