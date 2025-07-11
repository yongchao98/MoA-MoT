import sys

# It's good practice to increase the limit for int to string conversion
# for very large numbers, though not strictly necessary for the calculation itself.
sys.set_int_max_str_digits(0)

def find_special_prime():
    """
    Finds the smallest prime p in a given list such that 6^(p-1) is congruent to 1 modulo p^2.
    This is the condition for Z[p√6] NOT being the ring of integers of Q(p√6).
    """
    
    primes = [17, 383, 1093, 66161, 534851]
    
    print("Searching for the smallest prime p such that 6^(p-1) ≡ 1 (mod p^2)...")
    print("-" * 60)

    for p in primes:
        print(f"Checking prime p = {p}...")
        
        modulus = p * p
        
        # Using Python's built-in pow(base, exp, mod) is highly efficient for modular exponentiation.
        result = pow(6, p - 1, modulus)
        
        print(f"Calculating 6^({p-1}) mod {p}^2. Result = {result}")
        
        if result == 1:
            print(f"\n>>> Found the prime: p = {p}")
            print(f"For p = {p}, the condition 6^({p-1}) ≡ 1 (mod {p}^2) is met.")
            print(f"Therefore, for this prime, Z[√[{p}]{6}] is not the ring of integers of Q(√[{p}]{6}).")
            
            print("\nThe numbers in the final equation 6^({p-1}) ≡ 1 (mod {p}^2) are:")
            print(6)
            print(p - 1)
            print(1)
            print(p)
            print(2)
            
            # We found the smallest prime in the list, so we can stop.
            break
        else:
            print("Condition not met. Continuing search...\n")

find_special_prime()