import cypari2

def solve_elliptic_curve_l_function_problem():
    """
    This function solves the problem by computing the rank of the elliptic curve E,
    the leading coefficients of its twisted L-functions, and their sum.
    """
    # Initialize the PARI/GP environment via cypari2 with sufficient precision
    pari = cypari2.Pari()
    pari.default('realprecision', 50)

    # Define the elliptic curve E: y^2 + y = x^3 - x^2 - 10x - 20
    # using its standard Weierstrass coefficients [a1, a2, a3, a4, a6].
    # For the given equation, a1=0, a2=-1, a3=1, a4=-10, a6=-20.
    E = pari.ellinit([0, -1, 1, -10, -20])

    # Compute the rank 'r' of the Mordell-Weil group E(Q).
    # The ellrank function returns a vector; the rank is the first component.
    # The rank of this curve is 0.
    r_vector = pari.ellrank(E)
    r = int(r_vector[0])

    # Find the two primitive cubic Dirichlet characters of conductor 7.
    # A character is cubic if its order is 3. For conductor 7 (a prime),
    # all non-trivial characters are primitive.
    G = pari.znchargroup(7)
    all_chars = G[1]  # G[1] is the list of characters for the group
    cubic_chars = []
    for chi in all_chars:
        # The order of the character is the 4th element (index 3) in the vector
        if chi[3] == 3:
            cubic_chars.append(chi)
    
    chi1 = cubic_chars[0]
    chi2 = cubic_chars[1]

    # Compute the leading coefficients 'a' and 'b' of the twisted L-functions.
    # The order of the zero at s=1 for these twisted L-functions is 1.
    # Thus, the leading coefficient is the value of the first derivative L'(E, 1, chi).
    # We use lfun(E, chi, s, d) to compute the d-th derivative of L(E, s, chi) at s.
    a = pari.lfun(E, chi1, 1, 1)
    b = pari.lfun(E, chi2, 1, 1)

    # Calculate the total sum. Since chi2 is the conjugate of chi1, b is the
    # conjugate of a, and their sum r + a + b is a real number.
    total_sum = (r + a + b).real()
    
    # Print the equation with all the computed numbers as requested.
    # The complex numbers a and b are shown in their entirety.
    print(f"{r} + ({a}) + ({b}) = {total_sum}")

    # The final answer must be rounded to four decimal places.
    # We will provide this in the required format separately.

solve_elliptic_curve_l_function_problem()