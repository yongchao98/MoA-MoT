import math

def solve():
    """
    Solves the problem by calculating the sum based on the determined condition for pairs (i, j).
    
    The condition for a pair (i, j) to be in the set S is that the coprime pair
    (i', j') = (i/gcd(i,j), j/gcd(i,j)) must be one of (1,1), (1,2), or (2,1).

    The total sum is \sum_{(i,j) in S} 1/2^(i+j).
    We can rewrite any pair (i,j) as (d*i', d*j') where d = gcd(i,j).
    The sum can be grouped by d and (i', j'):
    Sum = \sum_{d=1 to inf} \sum_{(i',j') in {(1,1),(1,2),(2,1)}} 1/2^(d*(i'+j'))

    This breaks down into three separate sums over d for each allowed coprime pair.
    For (i',j')=(1,1): sum_d 1/2^(2d) = sum_d (1/4)^d
    For (i',j')=(1,2): sum_d 1/2^(3d) = sum_d (1/8)^d
    For (i',j')=(2,1): sum_d 1/2^(3d) = sum_d (1/8)^d

    The sum of a geometric series sum_{d=1 to inf} r^d is r/(1-r).
    """

    # For (i', j') = (1,1)
    r1 = 1/4
    sum1 = r1 / (1 - r1)
    print(f"For coprime pair (1,1), i'+j'=2. The geometric series sum is (1/4)/(1-1/4) = {sum1}")

    # For (i', j') = (1,2)
    r2 = 1/8
    sum2 = r2 / (1 - r2)
    print(f"For coprime pair (1,2), i'+j'=3. The geometric series sum is (1/8)/(1-1/8) = {sum2}")

    # For (i', j') = (2,1)
    r3 = 1/8
    sum3 = r3 / (1 - r3)
    print(f"For coprime pair (2,1), i'+j'=3. The geometric series sum is (1/8)/(1-1/8) = {sum3}")
    
    total_sum = sum1 + sum2 + sum3
    
    term1_num, term1_den = 1, 3
    term2_num, term2_den = 2, 7
    final_num = term1_num * term2_den + term2_num * term1_den
    final_den = term1_den * term2_den
    
    print(f"\nThe total sum is the sum of these three values:")
    print(f"{sum1:.4f} + {sum2:.4f} + {sum3:.4f} = {total_sum:.4f}")
    print(f"In fractions, this is {term1_num}/{term1_den} + {term2_num}/{term2_den} = {final_num}/{final_den}")
    print("\nThe final equation is:")
    print(f"(1/4)/(1-1/4) + (1/8)/(1-1/8) + (1/8)/(1-1/8) = 1/3 + 1/7 + 1/7 = 1/3 + 2/7 = (7+6)/21 = {final_num}/{final_den}")

solve()