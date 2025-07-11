import math

def solve():
    """
    This function calculates the satisfaction values for W_1 and W_2 and their ratio.
    """
    
    # Part 1: Calculate s(N, W1)
    # W1 is a core committee that minimizes satisfaction for group N.
    # It must contain all 8 candidates from C = {1,...,8}.
    # The remaining 12 members are chosen to minimize intersection with D_i = A(i)\C.
    # This means picking 8 candidates from {25,...,32} and 4 from {9,...,24}.
    
    # Satisfaction from the common part C
    s_W1_common = 8 * 8
    # Satisfaction from the unique parts D_i
    s_W1_unique = 4
    s_W1 = s_W1_common + s_W1_unique
    
    print("For committee W1 (core):")
    print(f"Satisfaction s(N, W1) = (8 voters * 8 common candidates) + 4 unique candidates")
    print(f"s(N, W1) = {s_W1_common} + {s_W1_unique} = {s_W1}")
    print("-" * 20)

    # Part 2: Calculate s(N, W2)
    # W2 is an EJR committee that minimizes satisfaction for group N.
    # The optimal construction to minimize s(N,W) has:
    # - 6 candidates from C = {1,...,8}
    # - 2 candidates from D_1 = {9,10} (to satisfy EJR for group N)
    # - 4 candidates from A(9) U A(10) (to satisfy EJR for voters 9 and 10)
    # - The remaining 8 slots are filled to minimize satisfaction for N.
    #   - 4 from the rest of A(9) U A(10)
    #   - 4 from the remaining D_i sets {D_2,...,D_8}

    # Satisfaction from the common part C
    s_W2_common_count = 6
    s_W2_common = 8 * s_W2_common_count
    
    # Satisfaction from the unique parts D_i
    # d_1 = 2 (to satisfy EJR for N), sum(d_i) for i>1 is 4 (from the filler candidates)
    s_W2_unique = 2 + 4
    
    s_W2 = s_W2_common + s_W2_unique
    
    print("For committee W2 (EJR):")
    print(f"Satisfaction s(N, W2) = (8 voters * {s_W2_common_count} common candidates) + {s_W2_unique} unique candidates")
    print(f"s(N, W2) = {s_W2_common} + {s_W2_unique} = {s_W2}")
    print("-" * 20)
    
    # Part 3: Calculate the ratio
    ratio = s_W1 / s_W2
    
    print("Final Ratio Calculation:")
    print(f"s(N, W1) / s(N, W2) = {s_W1} / {s_W2}")
    
    # The final answer in the required format
    print(f"\n<<<34/27>>>")
    

solve()