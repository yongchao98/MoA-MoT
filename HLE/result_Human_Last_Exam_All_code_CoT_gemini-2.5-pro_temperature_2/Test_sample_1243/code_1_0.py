import math

def solve():
    # Part (a)
    # The number of vectors is the sum of combinations C(n, k) for n=16 and k = 4, 8, 12, 16.
    n_a = 16
    c4_a = math.comb(n_a, 4)
    c8_a = math.comb(n_a, 8)
    c12_a = math.comb(n_a, 12)
    c16_a = math.comb(n_a, 16)
    # The final equation is 1820 + 12870 + 1820 + 1 = 16511
    val_a = c4_a + c8_a + c12_a + c16_a

    # Part (b)
    # This is a question of existence, which was demonstrated in the thinking steps.
    # An example vector is x = (4, 4, 2, 0, 0, 0, 0, 0), whose sum of squares is 36,
    # which is a multiple of 18. Its components are all even. It is 3-primitive.
    # Therefore, the resulting neighbor N_3(x) is even.
    val_b = "yes"

    # Part (c)
    # The number is the sum of combinations C(n, k) for n=12 and k = 4, 8, 12.
    n_c = 12
    c4_c = math.comb(n_c, 4)
    c8_c = math.comb(n_c, 8)
    c12_c = math.comb(n_c, 12)
    # The final equation is 495 + 495 + 1 = 991
    val_c = c4_c + c8_c + c12_c

    # Print the final answer in the requested format.
    # The breakdown of each calculation is included in the plan above.
    # For (a), the total is 1820 + 12870 + 1820 + 1 = 16511.
    # For (c), the total is 495 + 495 + 1 = 991.
    print(f"(a) {val_a}; (b) {val_b}; (c) {val_c}")

solve()