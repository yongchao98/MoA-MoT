import math
from collections import defaultdict

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
    factors = defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def main():
    # The proposed set of 16 distinct integers
    S = [2, 3, 4, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 30, 32, 36]
    N = max(S)

    print(f"The proposed set of 16 distinct numbers is:\n{S}")
    print(f"The maximum value N in this set is: {N}\n")

    # Verify the conditions for this set
    
    # Get all prime factors involved in the set S
    all_primes = set()
    for num in S:
        factors = get_prime_factorization(num)
        for p in factors:
            all_primes.add(p)

    print("Verifying conditions for each prime:")
    valid_set = True
    exponent_sums = {}

    for p in sorted(list(all_primes)):
        exponents = [get_prime_factorization(n)[p] for n in S]
        total_exponent_sum = sum(exponents)
        max_exponent = max(exponents)
        
        print(f"  For prime {p}:")
        print(f"    Exponents: {exponents}")
        print(f"    Sum of exponents: {total_exponent_sum}")

        # Check if sum is divisible by 4
        if total_exponent_sum % 4 != 0:
            print(f"    Condition FAILED: Sum of exponents ({total_exponent_sum}) is not divisible by 4.")
            valid_set = False
        else:
            print(f"    Condition PASSED: Sum of exponents ({total_exponent_sum}) is divisible by 4.")
            magic_sum = total_exponent_sum // 4
            exponent_sums[p] = magic_sum
            print(f"    Required magic sum for exponents (S_p): {magic_sum}")

            # Check if max exponent is not greater than the magic sum
            if max_exponent > magic_sum:
                print(f"    Condition FAILED: Max exponent ({max_exponent}) is greater than the magic sum ({magic_sum}).")
                valid_set = False
            else:
                print(f"    Condition PASSED: Max exponent ({max_exponent}) is not greater than the magic sum ({magic_sum}).")
    
    if valid_set:
        print("\nAll conditions are satisfied for the set with N=36.")
        
        # Calculate the magic product K
        K = 1
        for p, s_p in exponent_sums.items():
            K *= p**s_p
        
        print(f"\nThe magic product K for each row and column would be:")
        
        equation_parts = []
        for p, s_p in sorted(exponent_sums.items()):
            equation_parts.append(f"{p}^{s_p}")
        
        print(f"K = {' * '.join(equation_parts)} = {K}")

    else:
        print("\nThe set is not valid.")

if __name__ == "__main__":
    main()
