def solve_diagonal_harmonics():
    """
    Solves the three-part question about diagonal harmonic polynomials.
    """

    # Part a: Find the bi-degree of the terminal polynomial.
    # An sl(2) string starting with a lowest weight vector of bi-degree (a, b)
    # terminates with a highest weight vector of bi-degree (b, a).
    # Given starting bi-degree (a, b) = (4, 3).
    # The terminal bi-degree is (b, a) = (3, 4).
    a_val = 4
    b_val = 3
    answer_a = f"({b_val}, {a_val})"

    # Part b: Provide the condition for a polynomial to be a valid string starter.
    # A polynomial of bi-degree (a, b) can be a string starter (lowest weight vector)
    # only if a >= b.
    # We assume a model where the polynomial is constructed from Delta_n(X)
    # (bi-degree (n(n-1)/2, 0)) by applying b operators. Each operator,
    # indexed by r_i, adds 1 to the y-degree and subtracts r_i from the x-degree.
    # This results in a = n(n-1)/2 - sum(r_i) and y-degree = b.
    # The condition a >= b becomes: n(n-1)/2 - sum(r_i) >= b.
    # This can be written as: sum_{i=1 to b} r_i <= n(n-1)/2 - b.
    # We format this expression as a string, including the numbers.
    n_str = "n"
    n_minus_1_str = "n-1"
    two_str = "2"
    b_str_var = "b"
    one_str = "1"
    i_str = "i"
    r_i_str = "r_i"
    sum_part = f"sum_{{{i_str}={one_str}}}^{{{b_str_var}}} {r_i_str}"
    inequality_sign = "<="
    rhs_part = f"{n_str}({n_str}-{one_str})/{two_str} - {b_str_var}"
    answer_b = f"{sum_part} {inequality_sign} {rhs_part}"

    # Part c: Determine if a polynomial of bi-degree (5, 2) can be constructed
    # using operators E_{r,0} for r = 1, 2.
    # E_{r,0} is multiplication by p_r(Y). A polynomial P with y-degree 2 would be
    # a linear combination of terms like p_1(Y)^2 * Q_1(X) and p_2(Y) * Q_2(X).
    # Applying F = sum(x_i * d/dy_i) to such a P results in a non-zero polynomial
    # unless P itself is zero. A string starter must be non-zero.
    answer_c = "No"

    print(f"a) {answer_a}")
    print(f"b) {answer_b}")
    print(f"c) {answer_c}")

solve_diagonal_harmonics()