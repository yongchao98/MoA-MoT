import math

def simulate_chocolate_passing(chocolates):
    """
    Simulates one minute of the chocolate passing game.
    """
    n = len(chocolates)
    next_chocolates = [0] * n
    for i in range(n):
        # Person i receives from person i-1 (with wrap-around)
        prev_person_idx = (i - 1 + n) % n
        temp_sum = (chocolates[i] + chocolates[prev_person_idx]) / 2
        
        # If the sum is odd, add 1 to make it even.
        # The problem states if anyone *ends up* with an odd number, they take one.
        # This is equivalent to rounding an odd temporary sum up to the next even number.
        # Since all chocolates are even, temp_sum is an integer.
        if temp_sum % 2 != 0:
            next_chocolates[i] = int(temp_sum + 1)
        else:
            next_chocolates[i] = int(temp_sum)
            
    return next_chocolates

def analyze_statements(initial_chocolates, max_steps):
    """
    Analyzes the three statements based on a simulation.
    """
    n = len(initial_chocolates)
    history = {'c': [], 'h': [], 'l': [], 'd': []}

    # Initial state
    c = initial_chocolates
    print(f"Initial state c_0: {c}")
    print(f"n = {n}\n")

    # Run simulation
    for i in range(max_steps):
        history['c'].append(c)
        h = max(c)
        l = min(c)
        d = h - l
        history['h'].append(h)
        history['l'].append(l)
        history['d'].append(d)
        
        print(f"--- Analyzing State at i = {i} ---")
        print(f"c_{i}: {c}, h_{i}: {h}, l_{i}: {l}, d_{i}: {d}\n")
        
        c = simulate_chocolate_passing(c)

    # Verify statements for i=0
    i = 0
    d_i = history['d'][i]
    l_i = history['l'][i]
    
    print("--- Verifying Statements for i = 0 ---")
    
    # Statement 1
    print("\nStatement 1: For any i >= 0, d(i+m) < d(i) for some m < n.")
    s1_holds = False
    for m in range(1, n):
        d_im = history['d'][i+m]
        check = d_im < d_i
        print(f"- m={m}: d({i+m}) < d({i})?  {d_im} < {d_i} is {check}")
        if check:
            s1_holds = True
            break
    if not s1_holds:
        print("No m found for i=0. Statement 1 is false.")

    # Statement 2 (check for i > 0)
    print("\nStatement 2: After the i-th minute for some i > 0, l(i) < l(i-1).")
    s2_holds_anywhere = False
    for j in range(1, max_steps):
        l_j = history['l'][j]
        l_j_minus_1 = history['l'][j-1]
        check = l_j < l_j_minus_1
        # We only need to show one check to demonstrate our logic
        if j <= 3: # Print first few checks
             print(f"- i={j}: l({j}) < l({j-1})? {l_j} < {l_j_minus_1} is {check}")
        if check:
            s2_holds_anywhere = True
    if not s2_holds_anywhere:
        print("The condition l(i) < l(i-1) was never met. Statement 2 is false.")
    
    # Statement 3
    print("\nStatement 3: For any i >= 0, there exists some m < n such that l(i+m) > l(i).")
    s3_holds = False
    for m in range(1, n):
        l_im = history['l'][i+m]
        check = l_im > l_i
        print(f"- m={m}: l({i+m}) > l({i})? {l_im} > {l_i} is {check}")
        if check:
            s3_holds = True
            break
    if not s3_holds:
        print("No m found for i=0. Statement 3 is false.")

    print("\nConclusion: All three statements are false.")


if __name__ == '__main__':
    # Use an initial state that is already converged to test the "for any i" condition.
    # n must be > 1. Let's use n=4.
    n_people = 4
    # max_steps must be at least n to check the statements for i=0.
    initial_config = [10] * n_people
    analyze_statements(initial_config, max_steps=n_people)
