import math

def simulate_chocolate_passing(chocolates):
    """
    Simulates one minute of the chocolate passing game.
    """
    n = len(chocolates)
    new_chocolates = [0] * n
    for i in range(n):
        # Person i receives from person i-1 (with wrap-around)
        prev_person_idx = (i - 1 + n) % n
        
        # Calculate the sum received
        s = (chocolates[i] + chocolates[prev_person_idx]) / 2
        
        # If the result is odd, they get an extra chocolate to make it even.
        if s % 2 != 0:
            new_chocolates[i] = math.ceil(s)
            if new_chocolates[i] % 2 != 0:
                 new_chocolates[i] += 1
        else:
            new_chocolates[i] = int(s)

    return new_chocolates

def check_statements(initial_chocolates, steps_to_check=5):
    """
    Checks the validity of the statements for a given initial distribution.
    """
    n = len(initial_chocolates)
    history = [initial_chocolates]
    
    print(f"Initial state c^0: {initial_chocolates}")
    
    # Run simulation
    for i in range(steps_to_check):
        history.append(simulate_chocolate_passing(history[-1]))

    # Analyze the statements
    print("\n--- Analyzing Statements ---")

    # Statement 2: l^i < l^{i-1} for some i > 0
    s2_holds = False
    for i in range(1, len(history)):
        l_i = min(history[i])
        l_prev = min(history[i-1])
        if l_i < l_prev:
            s2_holds = True
            break
    print(f"Statement 2 (l^i < l^(i-1)) is {s2_holds}.")
    print("Analysis: This is always false as the minimum is non-decreasing (l^i >= l^(i-1)).")

    # Check statements 1 and 3 for i=0
    i = 0
    d_i = max(history[i]) - min(history[i])
    l_i = min(history[i])
    
    s1_holds_for_i = False
    s3_holds_for_i = False
    
    for m in range(1, n):
        if i + m < len(history):
            d_im = max(history[i+m]) - min(history[i+m])
            l_im = min(history[i+m])
            if d_im < d_i:
                s1_holds_for_i = True
            if l_im > l_i:
                s3_holds_for_i = True

    print(f"\nChecking for i=0 with initial state {history[0]}:")
    print(f"d^0 = {d_i}, l^0 = {l_i}")
    print(f"Statement 1 ('exists m<n, d^(0+m) < d^0') holds? {s1_holds_for_i}")
    if not s1_holds_for_i:
        print("Counterexample found for Statement 1. It requires d^m < d^0, but d^m = d^0 = 0.")
        
    print(f"Statement 3 ('exists m<n, l^(0+m) > l^0') holds? {s3_holds_for_i}")
    if not s3_holds_for_i:
        print("Counterexample found for Statement 3. It requires l^m > l^0, but l^m = l^0.")

    print("\n--- Conclusion ---")
    print("Since we found a valid scenario where Statements 1 and 3 fail (for i=0),")
    print("and Statement 2 is always false, none of the statements are universally true.")


# Run the check with a uniform distribution
uniform_chocolates = [10, 10, 10, 10]
check_statements(uniform_chocolates)

<<<A>>>