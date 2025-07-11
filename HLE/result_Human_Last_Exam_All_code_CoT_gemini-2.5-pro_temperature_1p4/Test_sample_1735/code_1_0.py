import math

def evaluate_permutation(p):
    """Calculates the cumulative sums and the max absolute cumulative sum for a permutation."""
    if not p:
        return 0, []
    
    current_sum = 0
    sums = []
    max_abs_sum = 0
    
    print("Evaluating sequence:", p)
    print("Cumulative sum calculation:")
    
    equation_parts = []
    for x in p:
        current_sum += x
        sums.append(current_sum)
        max_abs_sum = max(max_abs_sum, abs(current_sum))
        
        # Build and print the equation for the current step
        if not equation_parts:
            equation_parts.append(f"({x})")
        else:
            if x < 0:
                equation_parts.append(f" - {-x}")
            else:
                equation_parts.append(f" + {x}")
        
        print(f"|{''.join(equation_parts)}| = |{current_sum}| = {abs(current_sum)}")
        
    print(f"Maximum cumulative impact: {max_abs_sum}\n")
    return max_abs_sum, sums

def greedy_f_solver(a):
    """
    Implements the greedy strategy from statement F:
    At each step k, choose x_k to minimize |S_{k-1} + x_k|.
    """
    remaining = list(a)
    result = []
    current_sum = 0
    
    while remaining:
        best_next_num = None
        min_abs_sum = float('inf')
        
        for x in remaining:
            if abs(current_sum + x) < min_abs_sum:
                min_abs_sum = abs(current_sum + x)
                best_next_num = x
        
        result.append(best_next_num)
        current_sum += best_next_num
        remaining.remove(best_next_num)
        
    return result

def demonstrate_falsity_of_f():
    """
    Demonstrates that statement F is false using a counterexample.
    """
    print("Demonstrating that statement F is false.\n")
    
    a = [10, -11, -11, 13]
    print(f"Input sequence: {a}\n")

    # Run the greedy solver from statement F
    p_greedy = greedy_f_solver(a)
    m_greedy, _ = evaluate_permutation(p_greedy)
    
    # Show a better, manually found permutation
    p_better = [-11, 13, -11, 10]
    m_better, _ = evaluate_permutation(p_better)
    
    print("--- Conclusion ---")
    print(f"The greedy algorithm produced a sequence with max impact {m_greedy}.")
    print(f"A different sequence exists with a smaller max impact of {m_better}.")
    print("Therefore, the greedy strategy described in statement F is not optimal.")

if __name__ == '__main__':
    demonstrate_falsity_of_f()
