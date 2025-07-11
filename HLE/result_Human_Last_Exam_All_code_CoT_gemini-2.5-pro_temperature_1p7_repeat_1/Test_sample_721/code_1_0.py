import numpy as np
from scipy.special import comb
from scipy.integrate import quad

def calculate_sums_and_C(num_terms=30):
    """
    Calculates the approximate values for S_1c, S_2c, and |C|.
    f_0(z) = sum(c_n * z^n), c_n is non-zero for n=4k+1
    c_{4k+1} = (binom(2k, k)) / (4^k * (4k+1))
    """
    S1_c = 0.0
    S2_c = 0.0
    
    # Sum over k for terms n=4k+1
    for k in range(num_terms):
        n = 4 * k + 1
        
        # Calculate coefficient c_n
        # c_n = special.binom(2*k, k) / (4**k * n)
        # Using comb from scipy to avoid large numbers in factorial
        binom_val = comb(2 * k, k, exact=True)
        c_n = binom_val / (4**k * n)
        
        # Add to sums
        S1_c += c_n
        S2_c += n * c_n**2

    # Calculate |C| = integral from 0 to 1 of 2/sqrt(1+u^4) du
    integrand = lambda u: 2 / np.sqrt(1 + u**4)
    C_val, _ = quad(integrand, 0, 1)

    return S1_c, S2_c, C_val

def main():
    """
    Main function to perform the check for function 2.
    """
    S1_c, S2_c, C_val = calculate_sums_and_C()
    
    A_val = np.sqrt(2)

    lhs_final = A_val**2 * S2_c
    rhs_final = C_val + A_val * S1_c
    
    # We derived a simpler sufficient condition: (2-sqrt(2))*S1_c <= |C|
    lhs_sufficient = (2 - A_val) * S1_c
    rhs_sufficient = C_val

    print("Analysis for Function 2:")
    print("This function maps the unit disk to a square.")
    print("The problem reduces to checking the inequality for the specific coefficients.")
    print("\nLet f_0(z) = sum(c_n z^n) be the canonical map to a square.")
    print("Our function is f(z) = C + A * f_0(z) with |A|=sqrt(2).")
    print("The inequality sum(n*|a_n|^2) <= sum(|a_n|) becomes:")
    print("|A|^2 * sum(n*|c_n|^2) <= |C| + |A| * sum(|c_n|)\n")

    print(f"Numerical approximation with 30 terms in the series:")
    print(f"|A| = {A_val:.4f}")
    print(f"|C| = integral(...) approx {C_val:.4f}")
    print(f"sum(|c_n|) approx S1_c = {S1_c:.4f}")
    print(f"sum(n*|c_n|^2) approx S2_c = {S2_c:.4f}")
    
    print("\nChecking the inequality: |A|^2 * S2_c <= |C| + |A| * S1_c")
    print(f"LHS = {lhs_final:.4f}")
    print(f"RHS = {rhs_final:.4f}")
    print(f"Is LHS <= RHS? {lhs_final <= rhs_final}\n")

    print("Checking the simpler sufficient condition: (2 - |A|) * S1_c <= |C|")
    print(f"LHS = {lhs_sufficient:.4f}")
    print(f"RHS = {rhs_sufficient:.4f}")
    print(f"Is LHS <= RHS? {lhs_sufficient <= rhs_sufficient}\n")

    if lhs_final <= rhs_final:
        print("Conclusion for function 2: The inequality holds.")
    else:
        print("Conclusion for function 2: The inequality does not hold.")

if __name__ == "__main__":
    main()
