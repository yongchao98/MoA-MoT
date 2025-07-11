import sys

def f(n, m):
    """
    This function returns 1 if and only if the first player has a winning position
    with a probability strictly more than 50% in a 2D-Generalized NIM game on a
    randomly generated n x m binary matrix.

    The analysis shows this is true if and only if n + m > 2.
    """
    
    print(f"Evaluating for a game on a {n}x{m} matrix...")
    
    # The simplified condition derived from the probability analysis is n + m > 2.
    # We will now check this condition.
    
    # The final equation to check is: n + m > 2
    lhs_value = n + m
    rhs_value = 2
    
    print(f"The condition is: n + m > 2")
    print(f"Breaking down the equation:")
    print(f"n = {n}")
    print(f"m = {m}")
    print(f"n + m = {lhs_value}")
    print(f"The number to compare against is {rhs_value}.")
    
    # Perform the comparison
    is_winning_prob_gt_50 = lhs_value > rhs_value
    
    if is_winning_prob_gt_50:
        result = 1
        print(f"Since {lhs_value} > {rhs_value} is True, the probability is > 50%. The function returns {result}.")
    else:
        result = 0
        print(f"Since {lhs_value} > {rhs_value} is False, the probability is not > 50%. The function returns {result}.")
        
    return result

if __name__ == '__main__':
    # We assume n and m are positive integers for a matrix dimensions.
    # Example 1: n=1, m=1. Here n+m=2, so the condition n+m>2 is false.
    print("--- Example 1 ---")
    f(1, 1)
    
    print("\n--- Example 2 ---")
    # Example 2: n=2, m=3. Here n+m=5, so the condition n+m>2 is true.
    f(2, 3)

    print("\n--- Example 3 ---")
    # Example 3: n=5, m=5. Here n+m=10, so the condition n+m>2 is true.
    f(5, 5)
