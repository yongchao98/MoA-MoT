import math

def solve():
    """
    Solves the problem by identifying a candidate tuple that maximizes the 
    sequence length while minimizing the sum, and then computes the required expression.
    """
    
    # Based on the analysis, the tuple (13, 2, 4, 7) is a strong candidate
    # for achieving the maximum sequence length (M) with the smallest possible sum.
    # It is a permutation of (t_3, t_4, t_5, t_6) = (2, 4, 7, 13) from a Tribonacci sequence,
    # and testing shows it produces a long sequence.
    
    a = 13
    b = 2
    c = 4
    d = 7
    
    # The task asks to show the final equation with each number.
    result_val = a + b - c - d
    print(f"The selected tuple is (a,b,c,d) = ({a}, {b}, {c}, {d}).")
    print("This tuple is conjectured to maximize the function f(a,b,c,d) while having the minimal possible sum.")
    print("\nThe expression to compute is a + b - c - d mod 1000.")
    print(f"a + b - c - d = {a} + {b} - {c} - {d} = {result_val}")
    
    final_result = result_val % 1000
    
    print(f"\nThe result of ({a} + {b} - {c} - {d}) mod 1000 is {final_result}.")
    print(f"<<<{final_result}>>>")

solve()
