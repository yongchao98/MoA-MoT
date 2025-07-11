import math

def h_X_rank(k):
    """
    Computes the rank of H_{SO(4)}^k(X; Q), which is the coefficient of t^k in
    the Poincare series P_X(t) = (1+t^3) / (1-t^4)^2.
    P_X(t) = (1+t^3) * sum_{j=0 to inf} (j+1)t^{4j}
           = sum (j+1)t^{4j} + sum (j+1)t^{4j+3}
    """
    if k < 0:
        return 0
    if k % 4 == 0:
        j = k // 4
        return j + 1
    elif k % 4 == 3:
        j = (k - 3) // 4
        return j + 1
    else:
        return 0

def r_A_rank(k):
    """
    Computes the rank of the k-th equivariant cohomology group A^k.
    """
    if k < 0:
        return 0
    # Ranks of H_{SO(4)}^*(SO(4)) are h_G^0=1, h_G^1=2, h_G^2=1
    if k == 0:
        return 1
    if k == 1:
        return 2
    if k == 2:
        # r_2 = h_G^2 + h_X^0
        return 1 + h_X_rank(0)
    if k == 3:
        # r_3 = h_X^1
        return h_X_rank(1)
    if k >= 4:
        # r_k = h_X^{k-2}
        return h_X_rank(k - 2)
    return 0

def solve():
    """
    Calculates the total rank of A in degrees <= 100.
    """
    total_rank = 0
    
    # Low degree ranks (k=0, 1, 2, 3)
    sum_low_degrees = 0
    for k in range(4):
        sum_low_degrees += r_A_rank(k)

    # Higher degree ranks (k=4 to 100)
    sum_high_degrees = 0
    for k in range(4, 101):
        sum_high_degrees += r_A_rank(k)

    total_rank = sum_low_degrees + sum_high_degrees
    
    # We can also compute the sum for higher degrees analytically
    # k = 4m+1 => r_k = h_X(4m-1)=h_X(4(m-1)+3)=m. m=1..24
    # k = 4m+2 => r_k = h_X(4m)=m+1. m=1..24
    
    sum_4m_plus_1 = sum(m for m in range(1, 25)) # k=5..101. For k<=100, m=1..24.
    sum_4m_plus_2 = sum(m+1 for m in range(1, 25)) # k=6..102. For k<=100, m=1..24.
    
    analytical_sum_high = sum_4m_plus_1 + sum_4m_plus_2
    
    print("The total rank is the sum of ranks from degree 0 to 100.")
    print(f"Sum of ranks for degrees 0, 1, 2, 3 is: {r_A_rank(0)} + {r_A_rank(1)} + {r_A_rank(2)} + {r_A_rank(3)} = {sum_low_degrees}")
    print(f"Sum of ranks for degrees 4 to 100 is: {sum_high_degrees}")
    print(f"Total Rank = {sum_low_degrees} + {sum_high_degrees} = {total_rank}")

solve()