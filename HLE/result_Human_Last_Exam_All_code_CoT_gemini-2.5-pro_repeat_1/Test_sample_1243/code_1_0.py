import math

def solve():
    """
    Solves the three parts of the problem and prints the formatted answer.
    """

    # Part (a)
    # Count the number of non-zero binary vectors of length 16 with weight k divisible by 4.
    # The weights can be 4, 8, 12, 16.
    a_val = math.comb(16, 4) + math.comb(16, 8) + math.comb(16, 12) + math.comb(16, 16)
    
    # Part (b)
    # This is a theoretical question answered by finding an example. The answer is 'yes'.
    b_val = "yes"

    # Part (c)
    # Count the number of distinct unimodular 2-neighbors of Z^12.
    # For each vector v = x mod 2 with weight k in {4, 8}, there are 2 lattices.
    # For v with weight k=12, there is 1 lattice since 12 = 4 (mod 8).
    c_val = 2 * math.comb(12, 4) + 2 * math.comb(12, 8) + 1 * math.comb(12, 12)

    # Print the final answer in the required format.
    print(f"(a) {a_val}; (b) {b_val}; (c) {c_val}")
    print(f"\n# Calculation details for verification:")
    print(f"# Part (a): C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16) = {math.comb(16, 4)} + {math.comb(16, 8)} + {math.comb(16, 12)} + {math.comb(16, 16)} = {a_val}")
    print(f"# Part (b): An example vector x = (7, 7, 1, 3, 3, 3, 3, 3) in Z^8 can generate an even neighbor lattice. So the answer is yes.")
    print(f"# Part (c): 2*C(12, 4) + 2*C(12, 8) + C(12, 12) = 2*{math.comb(12, 4)} + 2*{math.comb(12, 8)} + {math.comb(12, 12)} = {c_val}")
    print("\n<<<" + f"(a) {a_val}; (b) {b_val}; (c) {c_val}" + ">>>")

solve()