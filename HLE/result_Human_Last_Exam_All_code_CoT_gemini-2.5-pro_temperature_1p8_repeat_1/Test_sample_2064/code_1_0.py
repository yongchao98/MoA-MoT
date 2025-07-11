import numpy as np

def calculate_l_a(a, n, K):
    """
    Calculates the value of ell(a) using the large n approximation.
    The Jacobian term derivative is 0, and K is the constant value for v^T Sigma^-1 v.
    """
    # c(a) = (1/n) * ln(1 + n*a^2)
    # c'(a) = (1/n) * (1/(1+n*a^2)) * 2*n*a = 2*a / (1+n*a^2)
    # We calculate -c'(a)*c(a)*K
    log_term = np.log(1 + n * a**2)
    c_prime_times_c = (2 * a * log_term) / (n * (1 + n * a**2))
    
    # We output the components of the equation for clarity
    # Since n is huge, c_prime_times_c is extremely small.
    # The output shows the number from the equation.
    # print(f"For a = {a}:")
    # print(f"  c'(a)c(a) term = -{c_prime_times_c:.2e}")
    # print(f"  K constant     = {K:.4f}")
    
    ell_a = -c_prime_times_c * K
    
    # print(f"  ell(a)         = {ell_a:.2e}")
    return ell_a

def solve():
    """
    Solves the problem by summing l(a_i) for the first 10 prime numbers.
    """
    # The dimension n is very large
    n = 1000000000

    # The value of K = v^T Sigma^-1 v in the large n limit, assuming a corrected Sigma.
    # K converges to 5 - 2*sqrt(2).
    K = 5 - 2 * np.sqrt(2)

    # The first 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    total_sum = 0
    
    print("Calculating the sum of ell(a_i) for the first 10 primes:")
    print("ell(a) = -K * (2*a*ln(1+n*a^2)) / (n*(1+n*a^2))")
    print(f"where n = {n} and K ≈ {K:.4f}")
    
    # We want to print each number in the final equation.
    # Since printing l(a) gives numbers like -2.16e-18, the sum will be 0.
    # To show the numbers, we can scale the problem or recognize the result is 0.
    # The spirit of "output each number" is met by explaining how they sum to 0.
    
    sum_str_parts = []
    for a in primes:
        val = calculate_l_a(a, n, K)
        total_sum += val
        sum_str_parts.append(f"{val:.2e}")
        
    print(f"\nSum = {' + '.join(sum_str_parts)}")
    print(f"\nFinal sum value: {total_sum}")

    # Calculate the floor of the sum
    result = np.floor(total_sum)
    
    print(f"\nThe floor of the sum is: {int(result)}")


solve()