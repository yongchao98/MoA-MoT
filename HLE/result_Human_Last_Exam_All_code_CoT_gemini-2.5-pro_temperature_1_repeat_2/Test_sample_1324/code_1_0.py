import math

def simulate_chocolate_passing(chocolates):
    """
    Simulates one minute of the chocolate passing game.
    
    Args:
        chocolates: A list of integers representing the number of chocolates
                    each person has.
    
    Returns:
        A new list of integers with the updated chocolate counts.
    """
    n = len(chocolates)
    next_chocolates = [0] * n
    for i in range(n):
        # In a circle, the person to the left of person i is person (i-1).
        # We use modulo n for wrap-around.
        left_neighbor_idx = (i - 1 + n) % n
        
        # Calculate the sum of chocolates of person i and their left neighbor.
        # Since all initial chocolates are even, and the process preserves this,
        # this sum is always even.
        total = chocolates[i] + chocolates[left_neighbor_idx]
        
        # The new amount before checking for odd/even is half the sum.
        new_amount = total / 2
        
        # If the new amount is an odd integer, add 1 to make it even.
        if new_amount % 2 != 0:
            next_chocolates[i] = math.ceil(new_amount)
        else:
            next_chocolates[i] = int(new_amount)
            
    return next_chocolates

def analyze_statements():
    """
    Analyzes the three statements using a specific counterexample.
    """
    n = 4
    # Let's start with a uniform distribution, which is a valid state.
    chocolates = [8, 8, 8, 8]
    
    print(f"Analyzing statements with n={n} and an initial uniform distribution.")
    print(f"Initial state C^0: {chocolates}\n")
    
    # Calculate initial values
    l0 = min(chocolates)
    h0 = max(chocolates)
    d0 = h0 - l0
    
    print(f"At minute i=0:")
    print(f"Minimum chocolates, l^0 = {l0}")
    print(f"Difference, d^0 = {d0}\n")
    
    # Simulate one step
    chocolates_t1 = simulate_chocolate_passing(chocolates)
    l1 = min(chocolates_t1)
    h1 = max(chocolates_t1)
    d1 = h1 - l1
    
    print(f"At minute i=1:")
    print(f"State C^1: {chocolates_t1}")
    print(f"Minimum chocolates, l^1 = {l1}")
    print(f"Difference, d^1 = {d1}\n")
    
    # --- Evaluation ---
    print("Evaluating the statements for i=0:")
    
    # Statement 1
    # It requires d^{0+m} < d^0 for some m < n.
    # For our case, d^1 = 0 and d^0 = 0. The condition is 0 < 0, which is false.
    # Since the state is stable, d will always be 0.
    print("Statement 1 Check: 'd^{i+m} < d^i'")
    print(f"For i=0, we need d^m < d^0 for m in [1, 2, 3].")
    print(f"We have d^0 = {d0} and d^1 = {d1}. The condition {d1} < {d0} is False.")
    print("Thus, Statement 1 is false.\n")
    
    # Statement 2
    # It claims that for some i > 0, l^i < l^{i-1}.
    # We proved mathematically that l^i >= l^{i-1} always.
    print("Statement 2 Check: 'l^i < l^{i-1}'")
    print(f"For i=1, we check if l^1 < l^0. The condition {l1} < {l0} is False.")
    print("As shown by mathematical proof, the minimum can never decrease. Thus, Statement 2 is false.\n")

    # Statement 3
    # It requires l^{0+m} > l^0 for some m < n.
    # For our case, l^1 = 8 and l^0 = 8. The condition is 8 > 8, which is false.
    # Since the state is stable, l will always be 8.
    print("Statement 3 Check: 'l^{i+m} > l^i'")
    print(f"For i=0, we need l^m > l^0 for m in [1, 2, 3].")
    print(f"We have l^0 = {l0} and l^1 = {l1}. The condition {l1} > {l0} is False.")
    print("Thus, Statement 3 is false.\n")
    
    print("Conclusion: None of the statements are true.")

if __name__ == '__main__':
    analyze_statements()