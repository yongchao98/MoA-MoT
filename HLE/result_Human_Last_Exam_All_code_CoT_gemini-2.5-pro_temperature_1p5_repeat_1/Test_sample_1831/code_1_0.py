import collections

def count_partitions(n):
    """
    Calculates the number of partitions of an integer n using dynamic programming.
    p(n) = p(n-1) + p(n-2) - p(n-5) - p(n-7) + ...
    """
    if n < 0:
        return 0
    dp = collections.defaultdict(int)
    dp[0] = 1
    for i in range(1, n + 1):
        k = 1
        while True:
            pentagonal1 = i - k * (3 * k - 1) // 2
            pentagonal2 = i - k * (3 * k + 1) // 2
            
            sign = (-1)**(k - 1)
            
            term1 = 0
            if pentagonal1 >= 0:
                term1 = sign * dp[pentagonal1]
            
            term2 = 0
            if pentagonal2 >= 0:
                term2 = sign * dp[pentagonal2]

            if term1 == 0 and term2 == 0:
                break
                
            dp[i] += term1 + term2
            k += 1
            
    return dp[n]

def solve():
    """
    Solves the problem by calculating the sum of partition numbers p(k) for k=1,2,3,4.
    """
    total_classes = 0
    for k in range(1, 5):
        pk = count_partitions(k)
        print(f"The number of partitions of {k}, p({k}), is {pk}.")
        total_classes += pk
    
    # Building the final equation string
    parts = [str(count_partitions(k)) for k in range(1, 5)]
    equation = " + ".join(parts)
    
    print(f"\nThe total number of equivalence classes represented is the sum:")
    print(f"{equation} = {total_classes}")

solve()