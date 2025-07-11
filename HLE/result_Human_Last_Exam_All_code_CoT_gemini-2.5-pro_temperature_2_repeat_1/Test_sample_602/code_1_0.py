import math

def solve_l(n: int) -> float:
    """
    Calculates the value of the function l(n) for n >= 5.

    The derivation shows that l(n) can be simplified to a direct formula involving n:
    l(n) = 2/n^2 * (n^2 + 1 - (2n-1)*sqrt(n^2 - n + 1))
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    # Calculate terms for the expression of l(n)
    n_squared = n * n
    sqrt_term = math.sqrt(n_squared - n + 1)
    
    term_in_parenthesis = n_squared + 1 - (2 * n - 1) * sqrt_term
    
    l_n = (2 / n_squared) * term_in_parenthesis
    
    return l_n

def print_equation(n: int):
    """
    Prints the components and the final value of the equation for l(n).
    """
    f1_A = 6
    a = math.sqrt(1 - (n - 1) / (n**2))
    b = 1/n

    S1 = (1/n**2) * (2*n**2 - 1 + (2*n-1)*math.sqrt(n**2 - n + 1))
    f1_M_Omega = 2 * S1

    l_n = solve_l(n)

    # Outputting numbers in the final equation as requested.
    print(f"For n = {n}:")
    print(f"The term f^(1)(A) evaluates to: {f1_A}")
    print(f"The term f^(1)(M*Omega) evaluates to: {f1_M_Omega}")
    print(f"The final equation is l({n}) = f^(1)(A) - f^(1)(M*Omega)")
    print(f"l({n}) = {f1_A} - {f1_M_Omega}")
    print(f"l({n}) = {l_n}")


# The user can execute this script. Since n is not specified, 
# we'll use n=5 as it's the smallest value in the function's domain.
n_val = 5
print_equation(n_val)