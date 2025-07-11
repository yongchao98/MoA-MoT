import math

def solve():
    """
    This function calculates the answers to the three parts of the problem.
    """
    # Part (a) calculation
    # Number of non-zero vectors in (Z/2Z)^16 with Hamming weight divisible by 4.
    # Weights can be k = 4, 8, 12, 16.
    c16_4 = math.comb(16, 4)
    c16_8 = math.comb(16, 8)
    c16_12 = math.comb(16, 12)
    c16_16 = math.comb(16, 16)
    ans_a_val = c16_4 + c16_8 + c16_12 + c16_16
    
    # Part (b) answer is based on logical deduction
    ans_b_val = "no"

    # Part (c) calculation
    # Number of non-zero vectors in (Z/2Z)^12 with Hamming weight divisible by 4.
    # Weights can be k = 4, 8, 12.
    c12_4 = math.comb(12, 4)
    c12_8 = math.comb(12, 8)
    c12_12 = math.comb(12, 12)
    ans_c_val = c12_4 + c12_8 + c12_12

    # Print the calculations as requested
    print(f"Calculation for (a): {c16_4} + {c16_8} + {c16_12} + {c16_16} = {ans_a_val}")
    print(f"Calculation for (c): {c12_4} + {c12_8} + {c12_12} = {ans_c_val}")
    
    # Print the final answer in the specified format
    print("\nFinal Answer:")
    final_answer = f"(a) {ans_a_val}; (b) {ans_b_val}; (c) {ans_c_val}"
    print(final_answer)

solve()