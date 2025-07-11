import functools
import operator

def get_grundy_value(a):
    """
    Calculates the Grundy value (g-number) for a single pile of size a,
    where a Zero-Move is available.
    """
    if a == 0:
        return 0
    if a % 2 == 1:  # a is odd
        return a + 1
    else:  # a is positive and even
        return a - 1

def solve_case(case_num, piles):
    """
    Solves a single case of the Nim game variant.
    """
    print(f"--- Case ({case_num}): a={piles} ---")
    
    grundy_values = [get_grundy_value(p) for p in piles]
    
    # Building the equation string
    g_expressions = []
    for p, g in zip(piles, grundy_values):
        g_expressions.append(f"g({p})={g}")
    print("Grundy values:", ", ".join(g_expressions))
    
    # Calculate Nim-sum
    nim_sum = functools.reduce(operator.xor, grundy_values)
    
    equation_parts = [str(g) for g in grundy_values]
    equation = " ^ ".join(equation_parts)
    print(f"Nim-sum = {equation} = {nim_sum}")
    
    if nim_sum == 0:
        winner = 'B'
        print("Result: Nim-sum is 0, Bob (second player) wins.\n")
    else:
        winner = 'A'
        print(f"Result: Nim-sum is {nim_sum} (non-zero), Alice (first player) wins.\n")
        
    return winner

def main():
    """
    Main function to solve all scenarios and print the final result.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]
    
    final_answer = ""
    for i, piles in enumerate(scenarios, 1):
        winner = solve_case(i, piles)
        final_answer += winner
        
    print(f"Final answer string: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    main()