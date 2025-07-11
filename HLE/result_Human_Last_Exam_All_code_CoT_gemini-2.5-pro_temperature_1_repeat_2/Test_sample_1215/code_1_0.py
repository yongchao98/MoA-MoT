def find_minimum_variables():
    """
    This function determines the minimum number of variables for a logically equivalent formula
    based on the given conditions.
    """
    print("Let's analyze the problem step by step.")
    print("Let phi be a formula with n variables (n >= 2).")
    print("1. phi is true for exactly 2^(n-1) out of 2^n total assignments.")
    print("2. We need to find the minimum number of variables, k, in a logically equivalent formula psi.")
    print("\nThe minimum number of variables in any equivalent formula is the number of 'essential' variables.")
    print("Let's find the minimum possible number of essential variables, k.")
    
    print("\nLet the boolean function for phi be f. If f has k essential variables, its truth value depends only on them.")
    print("Let T_k be the number of assignments to these k variables that make the function true.")
    print("The other (n-k) variables are non-essential ('dummy'), and their values can be chosen in 2^(n-k) ways.")
    print("The total number of satisfying assignments (out of 2^n) is N_true = T_k * 2^(n-k).")
    print("From the problem, we know N_true = 2^(n-1).")
    
    print("\nSo, we have the equation: T_k * 2^(n-k) = 2^(n-1)")
    print("Solving for T_k, we get: T_k = 2^(n-1) / 2^(n-k) = 2^((n-1) - (n-k))")
    print("The final equation is: T_k = 2^(k-1)")
    
    print("\nNow, we need to find the smallest non-negative integer k for which a boolean function with k variables")
    print("can have T_k = 2^(k-1) satisfying assignments.")
    print("A valid function must have 0 <= T_k <= 2^k, and T_k must be an integer.")
    
    k = 0
    while True:
        print(f"\n--- Checking k = {k} ---")
        
        # For k=0, the equation is T_0 = 2^(-1) = 0.5
        if k == 0:
            print("For k = 0, the equation is T_0 = 2^(0-1) = 2^(-1) = 0.5.")
            print("This is not an integer. A function cannot be true for half of an assignment.")
            print("Therefore, k cannot be 0.")
            k += 1
            continue

        # For k > 0, T_k = 2^(k-1) is an integer.
        t_k = 2**(k - 1)
        max_t_k = 2**k
        
        print(f"For k = {k}, the equation is T_{k} = 2^({k}-1) = {t_k}.")
        print(f"The total number of possible assignments for k variables is 2^k = {max_t_k}.")
        print(f"We must check if a function can exist, i.e., if 0 <= T_{k} <= 2^k.")
        print(f"Checking: 0 <= {t_k} <= {max_t_k}.")
        
        # This condition is always true for k >= 1.
        if 0 <= t_k <= max_t_k:
            print(f"The condition is met. A function with k={k} essential variables can exist.")
            print(f"For example, a function of one variable f(p_1) = p_1 is true for 1 assignment, which is 2^(1-1).")
            print(f"Such a function satisfies all the problem's conditions.")
            print(f"\nSince k=0 is not possible, the minimum possible value for k is {k}.")
            print("The fact that the formula psi can only use conjunction and negation does not change this,")
            print("because the set of connectives {AND, NOT} is functionally complete, meaning it can represent any boolean function,")
            print("including a simple function of one variable like p_1.")
            
            # Final answer format
            print(f"\n<<<1>>>")
            return
        
find_minimum_variables()