def solve():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of SO(4) \ X
    in degrees up to 100, where X is a 3-dimensional invariant submanifold.
    
    The Poincare series for the cohomology ring A is given by:
    P_A(t) = (t^3 + t^4 + t^6 + t^7) / (1 - t^4)^2
           = (t^3 + t^4 + t^6 + t^7) * sum_{k=0 to inf} (k+1)*t^(4k)
    
    The rank of A in degree n is the coefficient of t^n in this series.
    We need to sum these ranks for n from 0 to 100.
    """
    
    total_rank = 0
    
    # Contribution from t^3 * sum(...)
    # We need to sum (k+1) for all k where 4k+3 <= 100
    limit_k1 = (100 - 3) // 4
    sum1 = sum(k + 1 for k in range(limit_k1 + 1))
    print(f"Contribution from terms t^(4k+3):")
    print(f"  4k+3 <= 100  =>  k <= {limit_k1}")
    print(f"  Sum_{{k=0}}^{{{limit_k1}}} (k+1) = {sum1}")
    total_rank += sum1
    
    # Contribution from t^4 * sum(...)
    # We need to sum (k+1) for all k where 4k+4 <= 100
    limit_k2 = (100 - 4) // 4
    sum2 = sum(k + 1 for k in range(limit_k2 + 1))
    print(f"Contribution from terms t^(4k+4):")
    print(f"  4k+4 <= 100  =>  k <= {limit_k2}")
    print(f"  Sum_{{k=0}}^{{{limit_k2}}} (k+1) = {sum2}")
    total_rank += sum2
    
    # Contribution from t^6 * sum(...)
    # We need to sum (k+1) for all k where 4k+6 <= 100
    limit_k3 = (100 - 6) // 4
    sum3 = sum(k + 1 for k in range(limit_k3 + 1))
    print(f"Contribution from terms t^(4k+6):")
    print(f"  4k+6 <= 100  =>  k <= {limit_k3}")
    print(f"  Sum_{{k=0}}^{{{limit_k3}}} (k+1) = {sum3}")
    total_rank += sum3
    
    # Contribution from t^7 * sum(...)
    # We need to sum (k+1) for all k where 4k+7 <= 100
    limit_k4 = (100 - 7) // 4
    sum4 = sum(k + 1 for k in range(limit_k4 + 1))
    print(f"Contribution from terms t^(4k+7):")
    print(f"  4k+7 <= 100  =>  k <= {limit_k4}")
    print(f"  Sum_{{k=0}}^{{{limit_k4}}} (k+1) = {sum4}")
    total_rank += sum4
    
    print("\nFinal Equation:")
    print(f"Total Rank = {sum1} + {sum2} + {sum3} + {sum4} = {total_rank}")
    
    # The final answer is just the number
    # print(f"<<<{total_rank}>>>")

solve()