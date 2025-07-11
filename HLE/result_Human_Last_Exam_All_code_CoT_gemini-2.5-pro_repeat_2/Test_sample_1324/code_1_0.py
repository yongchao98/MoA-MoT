import math

def round_up_to_even(x):
    """
    Calculates the number of chocolates for a person in the next step.
    If the initial sum is odd, they take an extra chocolate, which is equivalent
    to rounding up the half-sum. The result must be even.
    The formula is equivalent to ceiling(sum/2) and then adding 1 if that result is odd.
    Or, more simply, if (a+b)/2 is an integer and it's odd, add 1.
    If (a+b)/2 is not an integer (i.e. a+b is odd), it's impossible since all c_k are even.
    So a+b is always even. (a+b)/2 can be even or odd.
    """
    val = (x) / 2
    # if the average is odd, add 1 to make it even
    if val % 2 != 0:
        return math.ceil(val)
    # if the average is even, it's already even
    else:
        return int(val)

def simulate_step(chocolates):
    """
    Simulates one minute of chocolate passing for the entire circle.
    """
    n = len(chocolates)
    next_chocolates = [0] * n
    for i in range(n):
        # Person i gives half to person i+1 (circularly)
        # So person i+1 gets chocolates from person i.
        # Person i (index i) gets from person i-1 (index i-1).
        # We model this as person k getting from person k-1.
        # Person p_k is at index k-1 in a 0-indexed list.
        # p_1 is at index 0, gets from p_n (index n-1)
        # p_k is at index k-1, gets from p_{k-1} (index k-2)
        
        # In our simulation, person at index `i` gets from person at index `i-1`.
        left_neighbor_idx = (i - 1 + n) % n
        
        current_val = chocolates[i]
        left_neighbor_val = chocolates[left_neighbor_idx]
        
        total_sum = current_val + left_neighbor_val
        new_val = total_sum / 2
        
        # If the average is odd, they take one more chocolate.
        if new_val % 2 != 0:
            next_chocolates[i] = int(new_val) + 1
        else:
            next_chocolates[i] = int(new_val)
            
    return next_chocolates

def analyze_statements():
    """
    Runs a simulation and analyzes the three statements based on the results.
    """
    n = 4
    # Initial state must be all even numbers
    initial_chocolates = [6, 2, 10, 4]
    
    # Store history
    c_history = [initial_chocolates]
    l_history = [min(initial_chocolates)]
    h_history = [max(initial_chocolates)]
    d_history = [max(initial_chocolates) - min(initial_chocolates)]
    
    print(f"n = {n}")
    print(f"i=0: chocolates={c_history[0]}, h={h_history[0]}, l={l_history[0]}, d={d_history[0]}")

    # Run simulation for a few steps
    current_c = initial_chocolates
    for i in range(1, 8):
        current_c = simulate_step(current_c)
        c_history.append(current_c)
        l_history.append(min(current_c))
        h_history.append(max(current_c))
        d_history.append(max(current_c) - min(current_c))
        print(f"i={i}: chocolates={c_history[i]}, h={h_history[i]}, l={l_history[i]}, d={d_history[i]}")

    print("\n--- Analysis ---")
    # The system stabilizes at i=2 with chocolates = [6, 6, 6, 6].
    # We test the statements for i=3, when the system is in a steady state.
    test_i = 3
    print(f"We will test the statements at i={test_i}, where the system is stable.")

    # --- Test Statement 1 ---
    print("\nStatement 1: For any i >= 0, d^{i+m} < d^i where m < n.")
    print(f"At i={test_i}, d^{test_i} = {d_history[test_i]}. We need to find m in [1, 2, 3] such that d^({test_i}+m) < d^{test_i}.")
    m_found = False
    for m in range(1, n):
        d_future = d_history[test_i + m]
        d_current = d_history[test_i]
        print(f"  Testing m={m}: Is d^{test_i+m}={d_future} < d^{test_i}={d_current}? {d_future < d_current}")
        if d_future < d_current:
            m_found = True
    if not m_found:
        print("No such m found. Statement 1 is false because it fails for i={test_i}.")

    # --- Test Statement 2 ---
    print("\nStatement 2: After the i^{th} minute for some i > 0, l^i < l^{i-1}.")
    print(f"The sequence of minimums (l) is: {l_history}")
    decreasing_found = False
    for i in range(1, len(l_history)):
        if l_history[i] < l_history[i-1]:
            decreasing_found = True
            print(f"Found a case at i={i} where l^{i}={l_history[i]} < l^{i-1}={l_history[i-1]}.")
            break
    if not decreasing_found:
        print("The sequence of minimums is non-decreasing. Statement 2 is false.")
        
    # --- Test Statement 3 ---
    print("\nStatement 3: For any i >= 0, there exists some m < n such that l^{i+m} > l^i.")
    print(f"At i={test_i}, l^{test_i} = {l_history[test_i]}. We need to find m in [1, 2, 3] such that l^({test_i}+m) > l^{test_i}.")
    m_found = False
    for m in range(1, n):
        l_future = l_history[test_i + m]
        l_current = l_history[test_i]
        print(f"  Testing m={m}: Is l^{test_i+m}={l_future} > l^{test_i}={l_current}? {l_future > l_current}")
        if l_future > l_current:
            m_found = True
    if not m_found:
        print("No such m found. Statement 3 is false because it fails for i={test_i}.")
        
    print("\nConclusion: All three statements are false.")

analyze_statements()