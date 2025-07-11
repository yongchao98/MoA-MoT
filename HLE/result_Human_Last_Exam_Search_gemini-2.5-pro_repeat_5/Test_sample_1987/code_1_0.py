memo = {}
def f_recursive(a):
    a = tuple(a)
    if a in memo:
        return memo[a]
    if a[0] < 0:
        return 0
    # Check if sorted
    for i in range(len(a) - 1):
        if a[i] > a[i+1]:
            memo[a] = 0
            return 0
    if sum(a) == 0:
        memo[a] = 1
        return 1
    
    res = 0
    for i in range(len(a)):
        b = list(a)
        b[i] -= 1
        res += f_recursive(tuple(b))
    
    memo[a] = res
    return res

# Calculate f(2, 4, 5)
result_2_4_5 = f_recursive((2, 4, 5))