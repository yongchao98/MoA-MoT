import sympy

def solve_hitting_time_moment():
    """
    Solves for the supremum of alpha such that E[tau^alpha] < infinity.
    The method involves calculating the mean drift of the Markov chain and applying a known theorem.
    """
    
    # Define c as a positive symbolic constant and n as the state variable
    c = sympy.symbols('c', positive=True)
    n = sympy.symbols('n')

    # Step 1: Define the transition probabilities for large n
    # p(n, n+1) is the probability of moving from n to n+1
    p_up = sympy.Rational(1, 2) - c / n
    # p(n, n-1) is the probability of moving from n to n-1
    p_down = sympy.Rational(1, 2) + c / n

    print("Step 1: The transition probabilities for large n are given as:")
    print(f"p(n, n+1) = {p_up}")
    print(f"p(n, n-1) = {p_down}")
    print("-" * 40)

    # Step 2: Calculate the mean drift mu_n = E[X_{k+1} - X_k | X_k = n]
    # The change in position is +1 with probability p_up and -1 with probability p_down.
    mu_n = (1 * p_up) + (-1 * p_down)
    
    # Simplify the expression for the drift
    mu_n_simplified = sympy.simplify(mu_n)

    print("Step 2: Calculate the mean drift 'mu_n'.")
    print(f"mu_n = (+1) * p(n, n+1) + (-1) * p(n, n-1)")
    print(f"mu_n = ({p_up}) - ({p_down})")
    print(f"The simplified mean drift is: mu_n = {mu_n_simplified}")
    print("-" * 40)

    # Step 3: Identify the constant 'a' from the asymptotic form mu_n ~ -a/n.
    # By comparing our calculated drift mu_n = -2*c/n with the form -a/n, we can find 'a'.
    a = -mu_n_simplified * n
    
    print("Step 3: A key theorem states that if the drift mu_n ~ -a/n,")
    print("then the supremum of alpha is a/2.")
    print("By comparing the calculated drift with the theoretical form, we find 'a'.")
    print(f"a = -n * mu_n = -n * ({mu_n_simplified})")
    print(f"So, a = {a}")
    print("-" * 40)
    
    # Step 4: Calculate the supremum of alpha, which is a/2.
    sup_alpha = a / 2
    
    print("Step 4: The supremum of alpha is calculated as a / 2.")
    print(f"sup{{alpha}} = {a} / 2")
    print(f"The final result is: {sup_alpha}")

if __name__ == "__main__":
    solve_hitting_time_moment()
