import sympy

def solve():
    """
    Calculates the 2-adic valuation of the discriminant of the polynomial
    f(x) = 2*x^5 + 2*x^3 + 1.
    """
    # Define the variable x and the polynomial f(x)
    x = sympy.Symbol('x')
    poly = 2*x**5 + 2*x**3 + 1

    # Calculate the discriminant of the polynomial
    d = sympy.discriminant(poly, x)

    print(f"The polynomial is f(x) = 2*x^5 + 2*x^3 + 1")
    print(f"The discriminant of f(x) is: {d}")

    # Calculate the 2-adic valuation of the discriminant
    # This is the number of times the discriminant is divisible by 2.
    if d == 0:
        # Should not happen for this polynomial
        valuation = float('inf')
        print("The 2-adic valuation is infinite.")
    else:
        n = int(d)
        valuation = 0
        if n % 2 == 0:
            temp_n = abs(n)
            while temp_n > 0 and temp_n % 2 == 0:
                valuation += 1
                temp_n //= 2
        print(f"To find the 2-adic valuation of {n}, we find the highest power of 2 that divides it.")
        print(f"{n} = 2^{valuation} * {n // (2**valuation)}")
        print(f"The 2-adic valuation v_2({n}) is {valuation}.")
        
    # The thickness is this valuation.
    thickness = valuation
    print(f"The thickness of the double point is equal to this valuation.")
    print(f"The final answer is {thickness}")


solve()