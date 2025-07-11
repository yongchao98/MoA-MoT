import math
from decimal import Decimal, getcontext

def solve_magnetization():
    """
    Calculates the minimum magnetization M_z(1) by searching over even n.
    The formula involves large numbers, so high-precision Decimal arithmetic is used.
    """
    # Set precision high enough to handle large intermediate values and cancellations
    getcontext().prec = 100

    min_mz = Decimal('inf')
    n_min = -1

    # Search for n_min over a range of even integers. A test run shows the minimum
    # occurs for n < 30.
    for n in range(2, 31, 2):
        n_dec = Decimal(n)
        pi_dec = Decimal(math.pi)

        # Calculate the constant part C_n = ((-2)^n) / (n! * pi^n * n^n)
        try:
            C_n_num = Decimal(-2)**n
            fact_n = Decimal(math.factorial(n))
            pi_n = pi_dec**n
            n_n = n_dec**n
            C_n_den = fact_n * pi_n * n_n
            C_n = C_n_num / C_n_den
        except (ValueError, OverflowError):
            continue

        # Calculate the sum part S
        # S = sum_{k=0}^{n+1} C(n+1, k) * P(4n, k) * (-1)^(n+1-k)
        total_sum = Decimal(0)
        for k in range(n + 2):
            try:
                comb_val = Decimal(math.comb(n + 1, k))

                # P(4n, k) = (4n)! / (4n-k)! can be very large.
                # Calculate it iteratively to maintain precision.
                perm_val = Decimal(1)
                for i in range(k):
                    perm_val *= Decimal(4 * n - i)

                term_sign = Decimal((-1)**(n + 1 - k))
                term = comb_val * perm_val * term_sign
                total_sum += term
            except (ValueError, OverflowError):
                total_sum = Decimal('inf') # Mark as invalid if overflow occurs
                break
        
        if total_sum == Decimal('inf'):
            continue

        current_mz = C_n * total_sum
        
        if current_mz < min_mz:
            min_mz = current_mz
            n_min = n
    
    # After finding n_min, print the detailed calculation for that n.
    print(f"The number of spins 'n' must be an even integer for the RHS to be positive.")
    print(f"Searching for the minimum M_z(1) over even n gives n_min = {n_min}.\n")
    print(f"The calculation for M_z(1) at n = {n_min} is based on the formula:")
    print(f"M_z(1) = [ ((-2)^n) / (n! * pi^n * n^n) ] * Sum")
    print(f"where Sum = sum_{{k=0 to n+1}} [ C(n+1, k) * P(4n, k) * (-1)^(n+1-k) ]\n")
    
    n = n_min
    n_dec = Decimal(n)
    pi_dec = Decimal(math.pi)
    
    # Recalculate and print the components for n_min
    C_n_num = Decimal(-2)**n
    fact_n = Decimal(math.factorial(n))
    pi_n = pi_dec**n
    n_n = n_dec**n
    C_n_den = fact_n * pi_n * n_n
    C_n = C_n_num / C_n_den
    
    total_sum = Decimal(0)
    for k in range(n + 2):
        comb_val = Decimal(math.comb(n + 1, k))
        perm_val = Decimal(1)
        for i in range(k):
            perm_val *= Decimal(4 * n - i)
        term_sign = Decimal((-1)**(n + 1 - k))
        total_sum += comb_val * perm_val * term_sign
        
    print(f"For n = {n_min}:")
    print(f"Constant Factor = ((-2)^{n}) / ({n}! * pi^{n} * {n}^{n}) = {C_n:.4e}")
    print(f"Sum Factor = {total_sum:.4e}")
    print(f"M_z(1) = {C_n:.4e} * {total_sum:.4e}")
    
    final_result = C_n * total_sum
    print(f"The minimum magnetization M_z(1) is {final_result:.8f}")
    
    print(f"\n<<<{final_result:.8f}>>>")

solve_magnetization()