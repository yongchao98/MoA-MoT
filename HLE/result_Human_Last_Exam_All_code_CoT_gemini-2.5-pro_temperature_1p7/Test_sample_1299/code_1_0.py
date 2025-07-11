import math

def get_prime_factorization(n):
    """Returns a dictionary of prime factors and their exponents."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def sigma_0(n_factors, power=1):
    """Calculates the number of divisors from a factorization dict."""
    if not n_factors:
        return 1
    res = 1
    expr_parts = []
    for p, e in n_factors.items():
        term = (power * e + 1)
        res *= term
        expr_parts.append(f"({power}*{e}+1)")
    
    if len(expr_parts) > 1:
      expr = "".join(expr_parts)
    else:
      expr = expr_parts[0]
      
    return res, expr

def calculate_expressions(ell):
    """Calculates and prints the expressions for |U_ell| and |T_ell|."""
    if not isinstance(ell, int) or ell <= 0:
        print("Please provide a positive integer for l.")
        return

    factors = get_prime_factorization(ell)
    exponents = list(factors.values())

    # Part A
    u_ell_val, u_ell_expr = sigma_0(factors, power=1)
    if ell == 1:
      print(f"A) |U_{ell}| = 1")
    else:
      print(f"A) |U_{ell}| = {u_ell_expr} = {u_ell_val}")

    # Part B
    if ell == 1:
        print("B) |T_1| = 1")
    else:
        t_ell_val, t_ell_expr = sigma_0(factors, power=2)
        print(f"B) |T_{ell}| = {t_ell_expr} - 1 = {t_ell_val} - 1 = {t_ell_val - 1}")

if __name__ == '__main__':
    # This block allows the user to test with a specific value for l.
    # The final output will be based on the logic within calculate_expressions.
    # As per the prompt, I will provide the expressions directly in the final answer.
    # To conform to the desired final output format, let's hardcode the solution expressions.
    
    # Let ell = p_1^{e_1} * ... * p_s^{e_s}
    # Part A expression
    part_a_expr = "A)(e_1+1)(e_2+1)...(e_s+1)"
    
    # Part B expression. My derivation showed |T_l| is sigma_0(l^2)-1 for l>1 and 1 for l=1.
    # This can be written as (2e_1+1)...(2e_s+1) - 1 for l>1.
    # The given `d` parameter does not arise from the definitions as stated.
    # The formula is the same for all l>1. The case l=1 is special.
    # There is no single expression using `e_i` that covers l=1 and l>1.
    # I'll provide the formula for l>1, as that is the general case involving e_i.
    part_b_expr = "B)((2*e_1)+1)((2*e_2)+1)...((2*e_s)+1) - 1, for l>1"
    
    # Final combined answer based on problem structure
    # A) |U_l| = sigma_0(l)
    # B) |T_l| = sigma_0(l^2) - 1 for l > 1 and |T_1|=1
    final_answer = "<<<A)sigma_0(l) B)(product_{i=1 to s}(2e_i+1)) - 1 + 2*delta_{l,1}>>>"
    
    # Since the final format requested is specific, let's re-format based on that.
    # Let's provide the answer in the most common symbolic forms.
    # For A: tau(l) or sigma_0(l) is standard.
    # For B: sigma_0(l^2) - 1 for l>1, and 1 for l=1.
    # A single expression that works for all l is sigma_0(l^2)-1+2*delta_{l,1} where delta is the Kronecker delta.
    # This seems to be the only way to reconcile the piecewise result with a single expression requirement.
    # The d parameter seems to be a red herring under the problem's explicit definitions.
    
    final_a = "sigma_0(l)"
    final_b = "(product_{i=1..s} (2*e_i + 1)) - 1 + 2*KroneckerDelta(l, 1)"
    # This is getting too complex. The simplest interpretation is that they are two separate questions.
    
    print("For a given l = p_1^e_1 * ... * p_s^e_s:")
    print("A) |U_l| has the expression: (e_1+1)(e_2+1)...(e_s+1)")
    print("B) |T_l| has the expression: ((2*e_1)+1)((2*e_2)+1)...((2*e_s)+1) - 1, for l > 1. For l=1, the value is 1.")
    print("\nThe problem is likely asking for these general symbolic forms.")
