import cypari2
import cmath

def solve_task():
    """
    This function solves the problem by computing the rank of the elliptic curve
    and the leading coefficients of its twisted L-functions to find r + a + b.
    """
    # Set up the PARI/GP environment with sufficient precision
    pari = cypari2.Pari()
    pari.set_real_precision(50)

    # 1. Define the elliptic curve E from its minimal Weierstrass equation
    # E: y^2 + y = x^3 - x^2 - 10x - 20 has coefficients [a1, a2, a3, a4, a6]
    # a1=0, a2=-1, a3=1, a4=-10, a6=-20
    E = pari.ellinit([0, -1, 1, -10, -20])

    # 2. Compute the rank r of the Mordell-Weil group E(Q)
    # The ellrank function returns a vector; the first component is the rank.
    r = pari.ellrank(E)[0]

    # 3. Define the two primitive cubic Dirichlet characters of conductor 7.
    # In PARI's Conrey notation, these are identified by [7, [2]] and [7, [4]].
    # These correspond to characters of order 3 in the character group of (Z/7Z)*, which is cyclic of order 6.
    chi1 = [7, [2]]
    chi2 = [7, [4]]

    # 4. Compute the leading coefficients a and b.
    # These are the first derivatives of the twisted L-functions at s=1.
    # We use the lfun(E, chi, s, d) function, where d is the derivative order.
    # The rank of these twisted L-functions is 1, so the leading term is the first derivative.
    a = pari.lfun(E, chi1, 1, 1)
    b = pari.lfun(E, chi2, 1, 1)

    # Convert PARI objects to standard Python types for printing
    r_val = int(r)
    a_val = complex(a)
    b_val = complex(b)

    # 5. Calculate the final sum r + a + b
    total_sum = r_val + a_val + b_val
    
    # As b is the complex conjugate of a, the sum is real.
    total_sum_real = total_sum.real
    
    # Round the final result to four decimal places
    final_answer = round(total_sum_real, 4)

    # Print each component of the sum and the final result.
    print(f"The rank of the elliptic curve is r = {r_val}")
    print(f"The leading coefficient 'a' for L(E, s, chi_1) at s=1 is: a = {a_val}")
    print(f"The leading coefficient 'b' for L(E, s, chi_2) at s=1 is: b = {b_val}")
    print(f"The required sum is r + a + b.")
    print(f"The equation with numerical values is: {r_val} + ({a_val}) + ({b_val})")
    print(f"The final numerical result is: {total_sum_real}")
    print(f"The result rounded to four decimal places is: {final_answer}")

if __name__ == "__main__":
    solve_task()