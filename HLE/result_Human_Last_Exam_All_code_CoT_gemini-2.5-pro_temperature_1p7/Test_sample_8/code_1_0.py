import sympy

def solve_conormal_regularity():
    """
    This function calculates the conormal regularity of R(sigma)f.
    
    The problem asks for the conormal space of R(sigma)f, where:
    - R(sigma) is the resolvent of the wave operator on the Schwarzschild metric.
    - f is a function in the conormal space A^(2+alpha)(X).
    """

    # We are given that f is in A^(2+alpha)(X).
    # The regularity index is s = 2 + alpha.
    alpha = sympy.Symbol('alpha')
    initial_regularity_s = 2 + alpha

    # The operator is R(sigma) = (Box_g - sigma^2)^-1.
    # Box_g is a second-order differential operator.
    # Therefore, (Box_g - sigma^2) is an elliptic operator of order m = 2.
    # The resolvent R(sigma) is its inverse, which is a pseudodifferential
    # operator of order k = -m.
    operator_order_k = -2

    # A pseudodifferential operator of order k maps the conormal space A^s to A^(s-k).
    # We calculate the new regularity index s'.
    final_regularity_s_prime = initial_regularity_s - operator_order_k

    # Print the explanation and the final result.
    print("Let the initial function f be in the conormal space A^s(X).")
    print(f"The initial regularity index is s = {initial_regularity_s}.")
    print("\nThe operator is the resolvent R(sigma), which is a pseudodifferential operator of order k.")
    print(f"The order of the operator is k = {operator_order_k}.")
    print("\nThe resulting function R(sigma)f belongs to the conormal space A^s'(X), where s' = s - k.")
    print("\nCalculating the new regularity index s':")
    print(f"s' = ({initial_regularity_s}) - ({operator_order_k})")
    
    # Show the calculation step by step
    calc_step1 = f"s' = 2 + {alpha} + 2"
    print(calc_step1)

    calc_step2 = f"s' = {2 + 2} + {alpha}"
    print(calc_step2)
    
    final_equation = f"s' = {final_regularity_s_prime}"
    print(final_equation)

    print("\nTherefore, R(sigma)f belongs to the conormal space:")
    print(f"A^({final_regularity_s_prime})(X)")

solve_conormal_regularity()
