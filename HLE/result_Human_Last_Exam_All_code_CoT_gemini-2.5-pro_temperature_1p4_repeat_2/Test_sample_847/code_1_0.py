import numpy as np

def solve():
    """
    This function calculates the constant C such that the minimal cost to sort the array
    is asymptotically C * n * ln(n).
    """

    # The optimal probability p_0 for a type 2 question is the real root of the
    # polynomial equation p^3 + p - 1 = 0.
    # We represent the polynomial by its coefficients: 1*p^3 + 0*p^2 + 1*p - 1.
    coeffs = [1, 0, 1, -1]

    # Find the roots of the polynomial.
    roots = np.roots(coeffs)

    # We are interested in the single real root, which lies between 0 and 1.
    p0 = 0.0
    for r in roots:
        if np.isreal(r):
            p0 = np.real(r)
            break

    # The minimal cost is asymptotically C * n * ln(n), where the constant C is -1/ln(p0).
    # This is the final equation for the constant C. The numbers involved are -1 and p0.
    constant_C = -1 / np.log(p0)
    
    # We output the numbers used in the final equation as per the instructions.
    # The final equation is: C = -1 / ln(p0)
    # The numbers are p0 and the calculated C.
    # print(f"The real root of p^3 + p - 1 = 0 is p0 = {p0:.5f}")
    # print(f"The constant C is -1 / ln({p0:.5f}) = {constant_C:.3f}")

    # The final answer is the value of the constant C.
    print(f"{constant_C:.3f}")

solve()