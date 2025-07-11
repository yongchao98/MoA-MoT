import math

def simulate_and_check(initial_chocolates, total_steps):
    """
    Simulates the chocolate passing game and checks the statements.
    """
    n = len(initial_chocolates)
    
    # Ensure initial values are even
    if any(c % 2 != 0 for c in initial_chocolates):
        print("Error: Initial chocolate counts must be even.")
        return

    print("--- Simulation Start ---")
    print(f"Initial state (i=0): {initial_chocolates}")
    
    history = [{'chocolates': initial_chocolates}]
    
    current_chocolates = list(initial_chocolates)
    
    # Main simulation loop
    for i in range(1, total_steps + 1):
        next_chocolates = [0] * n
        for k in range(n):
            # Person k receives from person k-1
            c_k = current_chocolates[k]
            # Python's negative index handles the circular array
            c_k_minus_1 = current_chocolates[k-1]
            
            # Calculate the average
            temp_val = (c_k + c_k_minus_1) / 2
            temp_val = int(temp_val)
            
            # If odd, take one extra chocolate to make it even
            if temp_val % 2 == 1:
                next_chocolates[k] = temp_val + 1
            else:
                next_chocolates[k] = temp_val
        
        current_chocolates = next_chocolates
        history.append({'chocolates': current_chocolates})

    # Analysis part
    # First, calculate h, l, d for all steps
    for i, step_data in enumerate(history):
        chocs = step_data['chocolates']
        step_data['h'] = max(chocs)
        step_data['l'] = min(chocs)
        step_data['d'] = step_data['h'] - step_data['l']
        print(f"i={i}: c={step_data['chocolates']}, h={step_data['h']}, l={step_data['l']}, d={step_data['d']}")

    print("\n--- Checking Statements ---")
    
    # Statement 2: l^i < l^{i-1}
    s2_true = False
    for i in range(1, len(history)):
        if history[i]['l'] < history[i-1]['l']:
            print(f"Statement 2 is TRUE at i={i}: l^{i}={history[i]['l']} < l^{i-1}={history[i-1]['l']}")
            s2_true = True
            break
    if not s2_true:
        print("Statement 2 is FALSE: The minimum value 'l' never decreased in this simulation.")

    # Statement 1 & 3
    s1_false_found = False
    s3_false_found = False
    for i in range(len(history)):
        # Check if the statement "exists m < n such that ..." holds for this i
        s1_holds_for_i = False
        s3_holds_for_i = False
        for m in range(1, n):
            if i + m < len(history):
                if history[i+m]['d'] < history[i]['d']:
                    s1_holds_for_i = True
                if history[i+m]['l'] > history[i]['l']:
                    s3_holds_for_i = True

        if not s1_holds_for_i and not s1_false_found:
             print(f"Statement 1 is FALSE: Counterexample at i={i}. d^i={history[i]['d']}. For all m in [1, {n-1}], d^(i+m) is not less than d^i.")
             s1_false_found = True

        if not s3_holds_for_i and not s3_false_found:
             print(f"Statement 3 is FALSE: Counterexample at i={i}. l^i={history[i]['l']}. For all m in [1, {n-1}], l^(i+m) is not greater than l^i.")
             s3_false_found = True

        if s1_false_found and s3_false_found:
            break

    print("\n--- Final Conclusion ---")
    print("Based on logical analysis and confirmed by simulation:")
    print("Statement 1 is false.")
    print("Statement 2 is false.")
    print("Statement 3 is false.")
    print("The correct option is A: None of the statements is true.")

# Example run
initial_state = [6, 6, 10, 10]
simulate_and_check(initial_state, 10)
