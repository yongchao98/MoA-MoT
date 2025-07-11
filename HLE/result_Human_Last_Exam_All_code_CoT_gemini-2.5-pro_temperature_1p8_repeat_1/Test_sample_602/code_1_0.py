import math

def solve_l_n(n):
    """
    Calculates the exact value of l(n) for a given integer n >= 5.

    The function follows the derived analytical solution.
    l(n) = 2 - c * (4 - 2/n) + 2/(n^2)
    where c = sqrt(n^2 - n + 1) / n.
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    # Calculate the value of c from the definition of the matrix M_n
    # c = sqrt(1 - (n-1)/n^2) = sqrt((n^2 - n + 1)/n^2)
    try:
        c_val = math.sqrt(n*n - n + 1) / n
    except ValueError:
        print("Error: Cannot compute square root of a negative number.")
        return

    # Derived analytical expression for l(n)
    # l(n) = 2 - 4*c + 2*c/n + 2/n^2
    
    term_4c = 4 * c_val
    term_2cn = 2 * c_val / n
    term_2n2 = 2 / (n * n)

    result = 2 - term_4c + term_2cn + term_2n2

    # Output the equation with numerical values, as requested
    print(f"For n = {n}:")
    print(f"The matrix M_n involves a constant c = sqrt(n^2 - n + 1) / n")
    print(f"c = sqrt({n*n} - {n} + 1) / {n} = sqrt({n*n - n + 1}) / {n} = {c_val:.10f}")
    print("\nThe function l(n) is calculated via the derived formula:")
    print("l(n) = 2 - 4*c + (2*c)/n + 2/(n^2)")
    print("\nSubstituting the values:")
    print(f"l({n}) = 2 - 4*({c_val:.10f}) + (2*{c_val:.10f})/{n} + 2/({n*n})")
    print(f"l({n}) = 2 - {term_4c:.10f} + {term_2cn:.10f} + {term_2n2:.10f}")
    print(f"l({n}) = {result:.10f}")


# Example calculation for n=5
n_value = 5
solve_l_n(n_value)
<<<l(n) = 2 - ( (4*n-2)*sqrt(n^2-n+1) )/n^2 + 2/n^2>>>