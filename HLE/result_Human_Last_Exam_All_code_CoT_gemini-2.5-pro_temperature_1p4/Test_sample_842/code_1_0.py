import sympy

def solve_and_check():
    """
    This function verifies the symmetry property f(t) = t^4 * f(t^-1)
    for each of the given answer choices.
    """
    t = sympy.Symbol('t')
    
    print("We derived the symmetry condition for f(t): f(t) = t^4 * f(1/t)")
    print("Now we test the given choices against this condition.")
    print("-" * 30)

    choices = {
        'A': 1,
        'B': t**2,
        'C': -1,
        'D': -t**3 + 3*t**2 - 2*t + 1,
        'E': 2*t**5 + 4*t**4 - 2*t**3 - 3*t**2
    }
    
    correct_answer = None

    for label, f_expr in choices.items():
        f_t = sympy.sympify(f_expr)
        f_t_inv = f_t.subs(t, 1/t)
        t4_f_t_inv = sympy.simplify(t**4 * f_t_inv)
        
        print(f"Testing Choice {label}:")
        # Print each number (coefficient) in the polynomial expression
        print(f"f(t) = {f_t}")
        # Show the calculation for the RHS of the symmetry condition
        print(f"t^4 * f(1/t) = t^4 * ({f_t_inv}) = {t4_f_t_inv}")
        
        if sympy.simplify(f_t - t4_f_t_inv) == 0:
            print(f"Result: Choice {label} satisfies the condition f(t) = t^4 * f(1/t).")
            correct_answer = label
        else:
            print(f"Result: Choice {label} does not satisfy the condition.")
        print("-" * 30)
        
    return correct_answer

if __name__ == '__main__':
    solve_and_check()