import sys

def solve_conormal_space():
    """
    Calculates the conormal space for the resolvent R(sigma) applied to a function f.
    
    This problem involves symbolic manipulation based on principles of microlocal analysis.
    - The input function f is in the conormal space A^{2+alpha}(X).
    - The operator R(sigma) is the resolvent of the wave operator Box_g.
    - Box_g is a second-order differential operator (order 2).
    - Its resolvent R(sigma) is a pseudo-differential operator of order -2.
    - A pseudo-differential operator of order 'm' maps a space A^s to A^{s-m}.
    """

    # Initial conormal order s = 2 + alpha
    initial_numeric_order = 2
    initial_symbolic_order = "alpha"

    # Order of the pseudo-differential operator R(sigma) is m = -2
    operator_order = -2

    # The final conormal order s' = s - m
    final_numeric_order = initial_numeric_order - operator_order
    
    print("The problem is to determine the conormal space of R(sigma)f.")
    print(f"The input function f belongs to A^s(X), where s = {initial_numeric_order} + {initial_symbolic_order}.")
    print("The resolvent R(sigma) is a pseudo-differential operator of order m = -2.")
    print("\nThe mapping rule for a pseudo-differential operator of order m is:")
    print("  Operator: A^s(X) -> A^(s-m)(X)")
    print("\nApplying this rule, the new conormal order is calculated as follows:")
    
    # Printing the equation with each number as requested
    print(f"  Final Order = ({initial_numeric_order} + {initial_symbolic_order}) - ({operator_order})")
    print(f"  Final Order = {initial_numeric_order} + {initial_symbolic_order} + {-operator_order}")
    print(f"  Final Order = {final_numeric_order} + {initial_symbolic_order}")
    
    final_space = f"A^{{{final_numeric_order}+{initial_symbolic_order}}}(X)"
    
    print(f"\nThus, R(sigma)f belongs to the conormal space: {final_space}")

    # The final answer must be wrapped in <<<>>>
    print(f"\n<<<{final_space}>>>", file=sys.stdout)

solve_conormal_space()