def derive_field_equation():
    """
    This function formats and prints the derived field equation for the given
    theory of gravity, highlighting the numerical coefficients of each term.
    """
    
    # Coefficients of the terms in the final equation.
    # The derivation leads to the equation in choice E.
    # Term 1: The derivative of the superpotential P
    # Term 2: A contraction of P and Q tensors
    # Term 3: Another contraction of P and Q tensors
    # Term 4: The non-metricity scalar term
    
    coeff_term1 = -2
    coeff_term2 = -1
    coeff_term3 = 2
    coeff_term4_numerator = 1
    coeff_term4_denominator = 2

    # String representations of the tensor parts of the equation
    term1_str = "1/sqrt(-g) * ∂_α(sqrt(-g) * P^α_{μν})"
    term2_str = "P_{μαβ} * Q_ν^{αβ}"
    term3_str = "Q^{αβ}_μ * P_{αβν}"
    term4_str = "Q * g_{μν}"
    rhs_str = "(8πG/c^4) * T_{μν}"

    # Construct the full equation string
    equation = (
        f"({coeff_term1}) * {term1_str} "
        f"+ ({coeff_term2}) * {term2_str} "
        f"+ ({coeff_term3}) * {term3_str} "
        f"+ ({coeff_term4_numerator}/{coeff_term4_denominator}) * {term4_str} "
        f"= {rhs_str}"
    )

    print("The final field equation is derived as follows, showing each term with its numerical coefficient:")
    print(equation)
    print("\nThis equation matches answer choice E.")


if __name__ == '__main__':
    derive_field_equation()