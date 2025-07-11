import math

def print_asymptotic_formula():
    """
    This function prints the derived asymptotic formula for d_B,delta.

    The problem asks for the asymptotic value of d_B,delta, the minimum degree of
    a polynomial p approximating e^(-x) within an error delta on the interval [0, B].
    This value is dependent on B and L = log(1/delta).

    The asymptotic behavior is different in two main regimes:
    1. For fixed B and delta -> 0 (L -> infinity), the theory of polynomial
       approximation for analytic functions shows that d = Theta(L / log(L)).
    2. For fixed delta and B -> infinity, the problem is significantly harder,
       as the polynomial must remain small over a large interval. Recent research
       shows that in this case, d = Theta(B^2).

    A single formula that captures the asymptotic behavior across the entire space
    where B + 1/delta -> infinity can be constructed by summing the dominant terms
    from these two regimes. The total degree must be sufficient to meet the
    requirements of both challenges.
    
    Thus, the combined asymptotic behavior is Theta(B^2 + L/log(L)).
    
    The code below will print this formula, explicitly stating the numerical
    coefficients as requested.
    """
    
    # Define the components of the formula as strings
    term1_coeff = "1"
    term1_base = "B"
    term1_exp = "2"
    
    term2_coeff = "1"
    term2_num = "L"
    term2_den = "log(L)"

    # Construct the final formula string
    # B**2 + L / log(L)
    final_formula = f"{term1_coeff} * {term1_base}**{term1_exp} + {term2_coeff} * {term2_num} / {term2_den}"
    
    print(final_formula)

print_asymptotic_formula()

# The final answer in the requested format
# Note: L is defined as log(delta^-1) in the problem.
print("<<<B**2 + L/log(L)>>>")