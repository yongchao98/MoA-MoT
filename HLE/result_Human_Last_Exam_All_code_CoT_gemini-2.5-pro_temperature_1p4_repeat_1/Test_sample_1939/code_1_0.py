import numpy as np
from scipy.special import comb

def F(a, b, n):
    """
    Calculates the number of integer points in aS_n + b(-S_n).
    Only implemented for n=2.
    """
    if n != 2:
        raise NotImplementedError("Only implemented for n=2 (d=3)")
    
    # Handle cases where combinations would be with negative upper indices
    # comb(n, k) = 0 if k > n or k < 0. scipy.special.comb handles this.
    term1 = comb(a, 2, exact=True)
    term2 = comb(b, 2, exact=True)
    
    return term1 + term2 + 2 * a * b + 2 * a + 2 * b + 1

def get_ehrhart_poly_coeffs(d):
    """
    Calculates the coefficients of the Ehrhart polynomial for a given dimension d.
    Only implemented for d=3.
    """
    if d != 3:
        raise NotImplementedError("Only implemented for d=3")
    
    n = d - 1
    # The degree of the polynomial is d. We need d+1 points to determine it.
    # We will calculate p(k) for k = 0, 1, 2, 3
    points = []
    for k in range(d + 1):
        pk = 0
        for j in range(k + 1):
            pk += F(k - j, j, n)
        points.append((k, pk))
    
    # Vandermonde matrix to solve for polynomial coefficients
    # c_d*k^d + ... + c_1*k + c_0 = p(k)
    x_vals = np.array([p[0] for p in points])
    y_vals = np.array([p[1] for p in points])
    
    V = np.vander(x_vals, d + 1)
    coeffs = np.linalg.solve(V, y_vals)
    return np.round(coeffs, 8)

def main():
    """
    Main function to calculate and display the roots of the Ehrhart polynomial for d=3.
    """
    d = 3
    # For d=3, we found p(z) = 1/3 * (z+1) * (2z^2+4z+3)
    # p(z) = 1/3 * (2z^3 + 4z^2 + 3z + 2z^2 + 4z + 3)
    # p(z) = 1/3 * (2z^3 + 6z^2 + 7z + 3)
    # Coeffs: [2/3, 2, 7/3, 1]
    
    # Let's derive the polynomial for p(z) and find its roots
    # p(z) = c_3*z^3 + c_2*z^2 + c_1*z + c_0
    
    # We found p(z) = 1/3 * (2z^3 + 6z^2 + 7z + 3). No, that expansion is wrong
    # (z+1)(2z^2+4z+3) = 2z^3 + 4z^2 + 3z + 2z^2 + 4z + 3 = 2z^3 + 6z^2 + 7z + 3. Oh, wait.
    # p(k) = 1/3 * (k+1)(2k^2+4k+3)
    # Let's check my summation again.
    # F(a,b,2) = 1/2 * (a^2-a+b^2-b + 4ab+4a+4b+2) = 1/2 * (a^2+b^2+4ab+3a+3b+2)
    # p(k) = sum over j=0..k of 1/2 * ( (k-j)^2+j^2+4(k-j)j + 3(k-j)+3j+2)
    # = 1/2 * sum( k^2-2kj+2j^2 + 4kj-4j^2 + 3k+2) = 1/2 * sum(k^2+3k+2 + 2kj-2j^2)
    # 1/2 * [ (k+1)(k^2+3k+2) + 2k*k(k+1)/2 - 2*k(k+1)(2k+1)/6 ]
    # = 1/2 * (k+1) [ (k+1)(k+2) + k^2 - k(2k+1)/3 ]
    # = 1/2 * (k+1) [ k^2+3k+2+k^2 - (2k^2+k)/3]
    # = (k+1)/6 * [3(2k^2+3k+2) - (2k^2+k)] = (k+1)/6 * [6k^2+9k+6-2k^2-k]
    # = (k+1)/6 * [4k^2+8k+6] = 1/3 * (k+1)(2k^2+4k+3). This is correct.
    
    coeffs = [2/3, 4/3, 3/3, 1/3] # p(z) = 1/3 * (2z^3+4z^2+3z+...) no..
    # 2z^3 + 6z^2 + 7z + 3
    coeffs_analytic = np.array([2/3, 6/3, 7/3, 3/3]) # [2/3, 2, 7/3, 1]
    
    roots = np.roots(coeffs_analytic)

    print("Analysis for d=3:")
    print("The Ehrhart polynomial is proportional to (z+1) * (2z^2 + 4z + 3).")
    print(f"The coefficients of the polynomial 2z^3 + 6z^2 + 7z + 3 are: {[2, 6, 7, 3]}")
    
    print("\nThe roots of the polynomial are:")
    for r in roots:
        print(f"{r:.6f}")
        
    print("\nThe real parts of the roots are:")
    for r in roots:
        print(f"{r.real:.6f}")

    print("\nFor d=2, the polynomial is (z+1)^2, with roots -1, -1. The real parts are -1.")
    print("For d=3, one root is -1 and the other two roots also have real part -1.")
    print("This suggests that for any d, every root of p has real part -1.")

if __name__ == '__main__':
    main()
