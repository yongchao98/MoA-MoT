import argparse

def pretty_print_formula(formula):
    """Prints a formula with reasonably placed parentheses."""
    # This is a simple heuristic for readability.
    # It wraps implications in parentheses if they are part of a tensor.
    if isinstance(formula, list):
        if formula[0] == 'tensor':
            return " ⊗ ".join(f"({pretty_print_formula(sub)})" if isinstance(sub, list) and sub[0] == 'limp' else pretty_print_formula(sub) for sub in formula[1:])
        elif formula[0] == 'limp':
            lhs = formula[1]
            rhs = formula[2]
            # Add parentheses around lhs if it's also an implication
            lhs_str = f"({pretty_print_formula(lhs)})" if isinstance(lhs, list) and lhs[0] == 'limp' else pretty_print_formula(lhs)
            return f"{lhs_str} ⊸ {pretty_print_formula(rhs)}"
    return formula

def generate_formulas(W, m, b):
    """
    Generates the linear logic formulas for the equipartitioning problem.
    
    Args:
        W (list[int]): The set of natural numbers.
        m (int): The number of partitions.
        b (int): The target sum for each partition.
    """
    
    # Check for the necessary condition
    if sum(W) != m * b:
        print(f"Error: The sum of elements in W ({sum(W)}) does not equal m * b ({m*b}).")
        print("An equipartition is impossible, so no logical encoding is generated.")
        return

    # 1. Define the state formulas S_0, S_1, ..., S_b
    # S_k are represented as nested lists for structure
    # ['limp', F1, F2] represents F1 ⊸ F2
    # ['tensor', F1, F2, ...] represents F1 ⊗ F2 ⊗ ...
    S = {0: '⊥'}
    for i in range(b):
        S[i+1] = ['limp', S[i], S[i]]

    print("--- State Formulas (S_k) ---")
    for i in range(b + 1):
        print(f"S_{i} = {pretty_print_formula(S[i])}")
    print("-" * 30)

    # 2. Define the function f(w)
    f_w_formulas = {}
    print("--- Resource Formulas (f(w)) ---")
    for w in sorted(list(set(W))):
        if w > b:
            print(f"Warning: w={w} is greater than b={b}. f({w}) is empty as it cannot be part of any partition.")
            f_w_formulas[w] = '1' # The identity for tensor
            continue
        
        links = []
        for i in range(b - w + 1):
            links.append(['limp', S[i], S[i+w]])
        
        if not links:
            f_w_formulas[w] = '1'
        elif len(links) == 1:
            f_w_formulas[w] = links[0]
        else:
            f_w_formulas[w] = ['tensor'] + links
        
        print(f"f({w}) = {pretty_print_formula(f_w_formulas[w])}")
    print("-" * 30)

    # 3. Define the goal formula C
    print("--- Goal Formula (C) ---")
    partition_goal = ['limp', S[0], S[b]]
    
    if m == 1:
        C = partition_goal
    else:
        C = ['tensor'] + [partition_goal] * m
    
    print(f"C = {pretty_print_formula(C)}")
    print("-" * 30)


if __name__ == '__main__':
    # Example usage: W = {1, 2, 3, 4}, m = 2, b = 5
    # The sum is 10, and m*b is 10. Partitions: {1,4} and {2,3}.
    example_W = [1, 2, 3, 4]
    example_m = 2
    example_b = 5

    print(f"Encoding EP(W, m, b) for W={example_W}, m={example_m}, b={example_b}\n")
    generate_formulas(example_W, example_m, example_b)
    
    final_answer_f = "f(w) = (S_0 ⊸ S_w) ⊗ (S_1 ⊸ S_{1+w}) ⊗ ... ⊗ (S_{b-w} ⊸ S_b), where S_0=⊥ and S_{k+1}=S_k ⊸ S_k"
    final_answer_C = "C = (S_0 ⊸ S_b) ⊗ ... ⊗ (S_0 ⊸ S_b) (m times)"
    
    # The final answer is the description of the function f and the formula C.
    # print(f"\n<<<f(w) = (S_0 ⊸ S_w) ⊗ (S_1 ⊸ S_{1+w}) ⊗ ... ⊗ (S_{b-w} ⊸ S_b) and C = (S_0 ⊸ S_b)^m, where S_0=⊥ and S_{k+1}=S_k ⊸ S_k>>>")
