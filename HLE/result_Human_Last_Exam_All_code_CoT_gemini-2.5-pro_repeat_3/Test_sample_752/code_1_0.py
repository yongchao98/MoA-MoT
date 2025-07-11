import collections

def number_of_partitions(n):
    """
    Calculates the number of partitions of an integer n using dynamic programming.
    """
    if n < 0:
        return 0
    if n == 0:
        return 1
    
    dp = collections.defaultdict(int)
    dp[0] = 1
    
    for i in range(1, n + 1):
        for j in range(i, n + 1):
            dp[j] += dp[j - i]
            
    return dp[n]

n = 9
result = number_of_partitions(n)
print(f"The number of partitions of {n} is p({n}).")
print(f"p(1) = 1")
print(f"p(2) = 2")
print(f"p(3) = 3")
print(f"p(4) = 5")
print(f"p(5) = 7")
print(f"p(6) = 11")
print(f"p(7) = 15")
print(f"p(8) = 22")
print(f"p(9) = {result}")
print(f"The number of conjugacy classes in the symmetric group S_{n} is given by the number of partitions of n.")
print(f"The cardinality of the coend is the number of conjugacy classes of Sigma_{n}.")
print(f"For n = {n}, the cardinality is p({n}) = {result}.")
