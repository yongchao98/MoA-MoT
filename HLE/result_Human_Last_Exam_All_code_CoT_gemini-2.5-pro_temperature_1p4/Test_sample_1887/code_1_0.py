import sys

def solve_set_theory_problem():
    """
    Solves the set theory problem by explaining the steps and printing the final answer.
    """
    
    # Use unicode for mathematical symbols if the terminal supports it
    try:
        sys.stdout.reconfigure(encoding='utf-8')
        omega = "\u03C9"
        aleph = "\u2135"
    except TypeError:
        omega = "omega"
        aleph = "aleph"

    print("### Step-by-Step Solution ###")
    print("\nStep 1: Characterize the set X of possible cofinalities.")
    print(f"Let kappa = 2^{omega}. We are given that kappa is a singular cardinal and kappa < {aleph}_({omega}_({omega}+5)).")
    print(f"The set X contains the possible values for lambda = cf(kappa).")
    print(f"  - By definition, lambda is a regular cardinal.")
    print(f"  - By König's Theorem, cf(2^{omega}) > {omega}, so lambda is an uncountable regular cardinal.")
    print(f"  - From the given constraints, cf(kappa) <= kappa < {aleph}_({omega}_({omega}+5)).")
    print(f"Therefore, X is the set of all uncountable regular cardinals strictly less than {aleph}_({omega}_({omega}+5)).")

    print("\nStep 2: Identify the uncountable regular cardinals.")
    print(f"A fundamental theorem in ZFC states that for any limit ordinal alpha > 0, the cardinal {aleph}_alpha is singular.")
    print(f"This implies that the only uncountable regular cardinals are successor cardinals, i.e., cardinals of the form {aleph}_(beta+1).")
    
    print("\nStep 3: Determine the order type of X.")
    print(f"The set X is {aleph}_(beta) where beta is a successor ordinal and beta < {omega}_({omega}+5).")
    print("The order type of X is the order type of its set of indices, which is the set of all successor ordinals less than "
          f"gamma = {omega}_({omega}+5).")
          
    print("\nStep 4: Conclude the order type.")
    print(f"For any limit ordinal gamma, the set of successor ordinals less than gamma has an order type of gamma.")
    print(f"The ordinal gamma = {omega}_({omega}+5) is a limit ordinal because its index, {omega} + 5, is greater than 0.")
    print(f"Thus, the order type of the set of indices is {omega}_({omega}+5).")

    # Define the numbers in the final expression
    number_5 = 5
    
    print("\n" + "="*40)
    print("Final Answer")
    print("="*40)
    print(f"The order type of X is: {omega}_({omega} + {number_5})")

solve_set_theory_problem()
