import math

def simulate_step(chocolates):
    """
    Simulates one minute of the chocolate passing game.
    
    Args:
        chocolates (list of int): The number of chocolates each person has.
    
    Returns:
        list of int: The new number of chocolates for each person after one minute.
    """
    n = len(chocolates)
    new_chocolates = [0] * n
    for k in range(n):
        # In the circle, the person to the left of p_k is p_{k-1}.
        # The index for the person to the left is (k-1) modulo n.
        k_minus_1 = (k - 1 + n) % n
        
        # Each person's new total is based on their own and their left neighbor's chocolates.
        # The sum is guaranteed to be even since all counts are even.
        total_from_pair = chocolates[k] + chocolates[k_minus_1]
        
        # The potential new amount is the average.
        prelim_amount = total_from_pair // 2
        
        # If the result is odd, they get one from the jar, making it even.
        if prelim_amount % 2 != 0:
            new_chocolates[k] = prelim_amount + 1
        else:
            new_chocolates[k] = prelim_amount
            
    return new_chocolates

def analyze_and_demonstrate():
    """
    Analyzes the three statements and uses a simulation to provide a counterexample for Statement 1.
    """
    print("### Step-by-step Analysis of the Statements ###")

    # --- Statement 2 Analysis ---
    print("\n--- Analysis of Statement 2 ---")
    print("Statement 2: 'After the i-th minute for some i > 0, l^i < l^(i-1)'.")
    print("This means the minimum number of chocolates, l^i, can decrease.")
    print("However, the new amount for any person k is c_k^(i+1) = ceil_even((c_k^i + c_(k-1)^i) / 2).")
    print("Since c_k^i >= l^i and c_(k-1)^i >= l^i, their average must be >= l^i.")
    print("The ceil_even function ensures c_k^(i+1) >= (c_k^i + c_(k-1)^i) / 2.")
    print("Therefore, c_k^(i+1) >= l^i for all k. This implies l^(i+1) >= l^i.")
    print("The minimum is non-decreasing. Thus, Statement 2 is FALSE.")

    # --- Statement 3 Analysis ---
    print("\n--- Analysis of Statement 3 ---")
    print("Statement 3: 'For any i >= 0, there exists some m < n such that l^(i+m) > l^i'.")
    print("This must hold for any initial state. Consider an equilibrium state where all c_k^0 are equal, e.g., c^0 = [8, 8, 8, 8].")
    print("In this case, l^0 = 8. At every subsequent step, c^i = [8, 8, 8, 8], so l^i = 8 for all i.")
    print("The condition l^(0+m) > l^0 becomes 8 > 8, which is false for all m.")
    print("Therefore, Statement 3 is FALSE.")

    # --- Statement 1 Analysis and Demonstration ---
    print("\n--- Analysis of Statement 1 ---")
    print("Statement 1: 'For any i >= 0, d^(i+m) < d^i where m < n'.")
    print("A strict interpretation implies this must hold for ALL m < n. We will show a counterexample for m=1.")
    print("Let's simulate with n=4 and an initial state c^0 = [10, 10, 10, 20].")
    
    chocolates_0 = [10, 10, 10, 20]
    n = len(chocolates_0)
    
    # Step i=1
    chocolates_1 = simulate_step(chocolates_0)
    h1 = max(chocolates_1)
    l1 = min(chocolates_1)
    d1 = h1 - l1
    print(f"\nAfter 1 minute (i=1):")
    print(f"c^1 = {chocolates_1}")
    print(f"d^1 = max(c^1) - min(c^1) = {h1} - {l1} = {d1}")

    # Step i=2
    chocolates_2 = simulate_step(chocolates_1)
    h2 = max(chocolates_2)
    l2 = min(chocolates_2)
    d2 = h2 - l2
    print(f"\nAfter 2 minutes (i=2):")
    print(f"c^2 = {chocolates_2}")
    print(f"d^2 = max(c^2) - min(c^2) = {h2} - {l2} = {d2}")

    print("\nNow we check Statement 1 for i=1 and m=1 (since m=1 < n=4).")
    print("The statement requires d^(1+1) < d^1, which is d^2 < d^1.")
    print("Substituting the values:")
    print(f"Is {d2} < {d1}?")
    print(f"The equation is {d2} < {d1}, which is {d2 < d1}.")
    print("Since the condition is not met, we have a counterexample.")
    print("Therefore, Statement 1 is FALSE.")
    
    print("\n### Final Conclusion ###")
    print("All three statements are false.")

# Execute the analysis and demonstration.
analyze_and_demonstrate()