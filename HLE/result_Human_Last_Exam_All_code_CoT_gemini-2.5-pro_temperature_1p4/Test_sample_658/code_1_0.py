import math

def analyze_asymptotic_degree(B, delta):
    """
    Analyzes the asymptotic value A(B, delta) for the minimum polynomial degree d_{B,delta}.

    The asymptotic relationship is d_{B,delta} = Theta(A(B, delta)), where
    A(B, delta) = B + L / log(L), with L = log(1/delta).

    This function calculates the value of A for given B and delta, and shows the
    contribution of each term to demonstrate the two main asymptotic regimes.
    """
    print(f"Analyzing for B = {B} and delta = {delta:.2e}")
    
    if delta <= 0 or delta >= 1:
        print("Error: delta must be in the interval (0, 1).")
        return
    if B < 1:
        print("Error: B must be >= 1.")
        return

    L = math.log(1 / delta)

    # The formula is asymptotic for large L. If L <= 1, log(L) is not positive.
    if L <= 1:
        print(f"Warning: L = log(1/delta) = {L:.4f} is not > 1. The asymptotic formula may not be reliable.")
        return
        
    term_B = B
    term_L_over_logL = L / math.log(L)
    
    total_asymptotic_value = term_B + term_L_over_logL

    print(f"  L = log(1/delta) = {L:.4f}")
    print(f"  The term 'B' contributes: {term_B:.4f}")
    print(f"  The term 'L/log(L)' contributes: {term_L_over_logL:.4f}")
    print(f"  The total asymptotic value A(B, delta) is: {total_asymptotic_value:.4f}\n")

def print_final_equation():
    """
    Prints the final symbolic equation and its coefficients.
    """
    print("The final asymptotic expression is of the form: c1 * B + c2 * L / log(c3 * L)")
    print("Based on the analysis, the coefficients are:")
    c1 = 1
    c2 = 1
    c3 = 1
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")


# --- Main Execution ---
# The script demonstrates the two asymptotic regimes with example values.

# Example 1: B is large relative to L/log(L). The 'B' term dominates.
print("Regime 1: B is the dominant term.")
analyze_asymptotic_degree(B=100, delta=1e-10)

# Example 2: L/log(L) is large relative to B. The 'L/log(L)' term dominates.
print("Regime 2: L/log(L) is the dominant term.")
analyze_asymptotic_degree(B=10, delta=1e-100)

# Finally, print the general formula and its coefficients as requested.
print_final_equation()