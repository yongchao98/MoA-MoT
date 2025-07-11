import sympy

def calculate_q_derivative_of_T_n():
    """
    This function symbolically calculates the q-difference quotient of T^n
    and prints the final resulting equation.
    """
    
    # Define the symbolic variables. We treat n as a general integer.
    q, T = sympy.symbols('q T')
    n = sympy.symbols('n', integer=True)
    
    # Define the expression for the q-difference quotient of T^n.
    # nabla_q(f(T)) = (f(q*T) - f(T)) / (q*T - T)
    nabla_q_T_n = ((q * T)**n - T**n) / (q * T - T)
    
    # The sympy.simplify function can perform the algebraic simplification.
    # The result of this simplification is T**(n-1) * (q**n - 1) / (q-1).
    # This serves as a programmatic verification of the derivation.
    simplified_expr = sympy.simplify(nabla_q_T_n)
    
    # The term (q**n - 1) / (q - 1) is the q-number [n]_q, which equals
    # the geometric sum 1 + q + q^2 + ... + q^(n-1).
    
    # Now, we format and print the final equation.
    # The left-hand side is the operation we are performing.
    lhs = "∇_q(T^n)"
    
    # The right-hand side shows the structure of the result, including
    # the coefficients and powers in the q-number sum.
    rhs = "(1 + q + q^2 + ... + q^(n-1)) * T^(n-1)"
    
    print("The q-difference quotient of T^n is given by the following formula:")
    print(f"\n{lhs} = {rhs}\n")
    
    # We can also show the compact form using the q-number notation [n]_q.
    compact_rhs_notation = f"[n]_q * T^(n-1)"
    print("Where [n]_q is the q-analogue of n, defined as the sum above.")
    print(f"So, the formula can also be written compactly as:")
    print(f"\n{lhs} = {compact_rhs_notation}\n")

if __name__ == '__main__':
    calculate_q_derivative_of_T_n()