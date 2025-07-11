import sympy as sp

def solve_and_print_formula():
    """
    This function derives the formula for P(n) by matching asymptotic series.
    """
    # Define the symbolic variables
    # L represents ln(n)
    # n is the variable
    L, n = sp.symbols('L n')

    # The Euler-Maclaurin formula gives the asymptotic expansion for ln(Q(n)).
    # We are interested in the remainder after the main terms (ln(A) + L**2/2).
    # log_Q_rem = L/(2*n) - (1-L)/(12*n**2) + (11-6*L)/(720*n**4) + O(L/n**6)
    # We will match our new formula's expansion against this.
    target_coeff_n2 = (L - 1) / 12
    target_coeff_n3 = 0  # The E-M expansion for this function has no n**-3 term

    # Define the structure of P(n) with unknown polynomial coefficients for L.
    P2 = sp.Symbol('P2') # Represents the numerator for the n**-2 term
    P3 = sp.Symbol('P3') # Represents the numerator for the n**-3 term
    P_n_symbolic = P2 / n**2 + P3 / n**3

    # Logarithm of the refinement part of the formula
    log_S_refinement = sp.log(1 + L/(2*n) + P_n_symbolic)

    # Expand as a series in 1/n. We need to go up to n**-3 to find P2 and P3.
    # The expansion is up to O(n**-4), so we ask for the series up to n=oo, limit=4
    series_log_S = log_S_refinement.series(n, sp.oo, 4).removeO()

    # --- Step 1: Solve for the n**-2 term (P2) ---
    # Extract the coefficient of n**-2 from the series
    coeff_n2 = series_log_S.coeff(n**-2)
    
    # Equate it to the target coefficient from the Euler-Maclaurin formula
    eq_n2 = sp.Eq(coeff_n2, target_coeff_n2)
    
    # Solve for P2
    P2_solved = sp.solve(eq_n2, P2)[0]

    # --- Step 2: Solve for the n**-3 term (P3) ---
    # Substitute the solution for P2 back into the series
    series_with_P2_solved = series_log_S.subs(P2, P2_solved)
    
    # Extract the coefficient of n**-3
    coeff_n3 = series_with_P2_solved.coeff(n**-3)
    
    # Equate it to the target coefficient (which is 0)
    eq_n3 = sp.Eq(coeff_n3, target_coeff_n3)
    
    # Solve for P3
    P3_solved = sp.solve(eq_n3, P3)[0]

    # --- Step 3: Construct and print the final formula for P(n) ---
    P_n_final = (P2_solved / n**2) + (P3_solved / n**3)
    
    print("The formula for P(n) is:")
    # Using pretty print for a clear mathematical representation
    print(sp.pretty(P_n_final, use_unicode=False, wrap_line=False))

# Run the function
solve_and_print_formula()