def solve_controller():
    """
    This function prints the derived set of all proper stabilizing controllers H_2(s).
    The controller is parametrized by a function K(s), which can be any stable
    and proper rational function.
    """

    # Numerator of H_2(s) is N_c(s) = (14s^2 + 41s + 26) + (s^2 - 1)K(s)
    # The numbers in the polynomial part are 14, 41, 26.
    # The numbers in the K(s) part are 1, -1.
    num_poly_part = "({}*s**{} + {}*s + {})".format(14, 2, 41, 26)
    num_k_part = "K(s)*({}*s**{} - {})".format(1, 2, 1)
    numerator = f"{num_poly_part} + {num_k_part}"

    # Denominator of H_2(s) is D_c(s) = (s^2 - 6s - 16) - sK(s)
    # The numbers in the polynomial part are 1, -6, -16.
    # The number in the K(s) part is -1.
    den_poly_part = "({}*s**{} - {}*s - {})".format(1, 2, 6, 16)
    den_k_part = "K(s)*s"
    denominator = f"{den_poly_part} - {den_k_part}"

    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print("H_2(s) = N(s)/D(s), where K(s) is any stable and proper rational function.\n")
    print("N(s) = {}".format(numerator))
    print("D(s) = {}".format(denominator))

solve_controller()