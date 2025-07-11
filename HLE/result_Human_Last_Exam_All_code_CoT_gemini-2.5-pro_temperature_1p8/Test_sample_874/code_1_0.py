import math

def solve():
    """
    Finds the specific tuple (a, b, c, d) and computes the required expression.
    """
    # The problem of finding the maximum length of a Ducci sequence is a known difficult problem.
    # The record-holders for length are constructed from Tribonacci-like sequences.
    # The general form is k*(0, T_n, T_{n+1}, T_{n+2}), where T_i are Tribonacci numbers.
    # To maximize the length, we must maximize n. To minimize the sum, k must be 1.

    # 1. Generate Tribonacci numbers (T_0=0, T_1=0, T_2=1) until we exceed the limit.
    limit = 10_000_000
    trib = [0, 0, 1]
    while trib[-1] <= limit:
        next_trib = trib[-1] + trib[-2] + trib[-3]
        trib.append(next_trib)

    # The last number is over the limit. The sequence generating tuple will use the three previous ones.
    # The maximal n for T_{n+2} <= limit is n=27.
    # The terms are T_27, T_28, T_29.
    
    t_27 = trib[27]
    t_28 = trib[28]
    t_29 = trib[29]

    # 2. Choose a canonical tuple (lexicographically first) that generates a max-length sequence.
    # This corresponds to (a, b, c, d) = (0, T_27, T_28, T_29)
    a = 0
    b = t_27
    c = t_28
    d = t_29

    # 3. Perform the final calculation: (a + b - c - d) mod 1000
    expression_val = a + b - c - d
    result = expression_val % 1000

    # 4. Print the explanation and result as requested.
    print("The tuple (a, b, c, d) that maximizes the Ducci sequence length under the given constraints")
    print("is based on the Tribonacci sequence. To minimize the sum, we choose the base sequence.")
    print("The lexicographically smallest optimal tuple is chosen:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print(f"d = {d}")
    print("\nCalculating the expression (a + b - c - d) mod 1000:")
    print(f"({a} + {b} - {c} - {d}) mod 1000")
    print(f"= {expression_val} mod 1000")
    print(f"= {result}")

solve()
print("<<<923>>>")