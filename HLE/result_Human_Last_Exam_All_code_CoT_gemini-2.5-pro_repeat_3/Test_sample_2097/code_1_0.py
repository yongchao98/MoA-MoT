import math
from scipy.special import comb

def p(n, k):
    """
    Calculates the permutation P(n, k) = n! / (n-k)! using integers
    to avoid floating point issues.
    """
    if k < 0 or k > n:
        return 0
    if k == 0:
        return 1
    res = 1
    for i in range(k):
        res *= (n - i)
    return res

def calculate_mz1(n):
    """
    Calculates M_z(1) for a given number of spins n.
    """
    # Calculate f(1)
    f1 = 0
    for j in range(n + 1):
        term = comb(n, j, exact=True) * ((-1)**j) * p(4 * n, j)
        f1 += term
        
    # Calculate f'(1)
    f1_prime = 0
    for j in range(n + 1):
        term = comb(n, j, exact=True) * ((-1)**j) * p(4 * n, j + 1)
        f1_prime += term
        
    diff = f1_prime - f1
    
    # Calculate the constant C_n
    try:
        cn = (2**n * n**(-n)) / (math.factorial(n) * math.pi**n)
    except (ValueError, OverflowError):
        # Use log for large n to maintain precision
        log_cn = n * math.log(2) - n * math.log(n) - math.lgamma(n + 1) - n * math.log(math.pi)
        cn = math.exp(log_cn)

    mz1 = cn * diff
    return mz1, f1, f1_prime

def solve():
    """
    Finds the minimum magnetization M_z(1) by searching over n.
    """
    min_mz1 = float('inf')
    n_min = 0
    
    # Search for n_min in a reasonable range, e.g., 1 to 15.
    # The magnitude of M_z(1,n) grows exponentially, so the minimum will be at small n.
    for n in range(1, 16):
        mz1, _, _ = calculate_mz1(n)
        if mz1 < min_mz1:
            min_mz1 = mz1
            n_min = n
            
    # Recalculate values for n_min to present them
    final_mz1, f1_val, f1_prime_val = calculate_mz1(n_min)
    
    diff_val = f1_prime_val - f1_val
    
    # Numerator and denominator of the constant C_n for presentation
    cn_num = 2**n_min
    cn_den_n_pow = n_min**n_min
    cn_den_fact = math.factorial(n_min)
    
    print(f"The minimum magnetization M_z(1) is found for n = {n_min} spins.")
    print("The final equation for M_z(1) is:")
    print(f"M_z(1) = ({cn_num} / ({cn_den_fact} * {cn_den_n_pow} * pi^{n_min})) * ({diff_val})")
    print("\nCalculated value:")
    print(f"M_z(1) = {final_mz1}")

solve()
<<< -11.69835011033282 >>>