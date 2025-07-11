import sympy

def find_infimum_c():
    """
    This function calculates the infimum of c for which the given Markov chain is transient.
    """
    # Step 1: Define symbolic variables for k (state) and c (parameter).
    k, c = sympy.symbols('k c')

    # Step 2: Define the transition probabilities and the corresponding jump sizes for large k.
    # P_k,j is the probability of moving from state k to state j.
    # The jump is j - k.
    transitions = {
        -2: sympy.Rational(1, 4),            # P_{k,k-2}
        -1: sympy.Rational(1, 4) - c/k,    # P_{k,k-1}
         1: sympy.Rational(1, 4) + c/k,    # P_{k,k+1}
         2: sympy.Rational(1, 4)            # P_{k,k+2}
    }

    # Step 3: Calculate the mean drift mu_k = E[X_{n+1} - k | X_n = k].
    # This is the expected jump size from state k.
    mu_k = sum(jump * prob for jump, prob in transitions.items())
    mu_k = sympy.simplify(mu_k)
    
    # The drift is of the form A/k. Find A = lim_{k->inf} (k * mu_k).
    A = sympy.limit(k * mu_k, k, sympy.oo)

    # Step 4: Calculate the second moment of the jump, S_k^2 = E[(X_{n+1} - k)^2 | X_n = k].
    s_k_sq = sum(jump**2 * prob for jump, prob in transitions.items())
    s_k_sq = sympy.simplify(s_k_sq)
    
    # The second moment tends to a constant B as k -> inf. Find B = lim_{k->inf} S_k^2.
    B = sympy.limit(s_k_sq, k, sympy.oo)
    
    print("Analysis of the Markov Chain for large k:")
    print(f"1. The mean drift (expected jump) from state k is: mu_k = {mu_k}")
    print(f"   For large k, this behaves like A/k, where A = lim(k * mu_k) = {A}")
    print(f"2. The second moment of the jump from state k is: S_k^2 = {s_k_sq}")
    print(f"   For large k, this tends to a constant B = lim(S_k^2) = {B}")
    
    # Step 5: Apply the Lamperti-Pakes criterion for transience.
    # The chain is transient if A > B/2.
    # The infimum of c is found at the boundary case, A = B/2.
    print("\nThe criterion for transience is A > B/2.")
    print("The infimum of c is found by solving the boundary equation A = B/2.")
    
    # Step 6: Solve the equation A = B/2 for c and print the steps.
    print("\nSolving the equation:")
    # Using floating point numbers for clear output of the calculation.
    A_val_expr = A
    B_val_num = float(B)
    
    print(f"   Equation: {A_val_expr} = {B_val_num} / 2")
    
    rhs_val = B_val_num / 2
    print(f"   => {A_val_expr} = {rhs_val}")
    
    # The expression for A is 2*c
    # We solve 2*c = rhs_val
    critical_c = rhs_val / 2
    print(f"   => c = {rhs_val} / 2")
    print(f"   => c = {critical_c}")
    
    print(f"\nThus, the chain is transient for c > {critical_c}.")
    print(f"The infimum of the set of such c is {critical_c}.")

# Execute the function to find the answer.
find_infimum_c()

print("\n<<<0.625>>>")