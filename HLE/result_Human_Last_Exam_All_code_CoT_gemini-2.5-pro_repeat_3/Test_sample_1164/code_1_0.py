def solve():
    """
    Finds and verifies the smallest integer n >= 2 that satisfies the problem's conditions.
    """
    p2_9 = 2**9
    p2_10 = 2**10
    p5_9 = 5**9
    p5_10 = 5**10

    n = 2
    while True:
        # Condition 1: The last 9 digits of n^k are eventually constant.
        # This is true if n^k(n-1) is divisible by 10^9 = 2^9 * 5^9 for large k.
        cond1_p2 = (n % 2 == 0) or ((n - 1) % p2_9 == 0)
        cond1_p5 = (n % 5 == 0) or ((n - 1) % p5_9 == 0)
        condition1_met = cond1_p2 and cond1_p5

        if condition1_met:
            # Condition 2: The last 10 digits of n^k are NOT eventually constant.
            # This is true if the condition for 10^10 stability is false.
            stable_10_p2 = (n % 2 == 0) or ((n - 1) % p2_10 == 0)
            stable_10_p5 = (n % 5 == 0) or ((n - 1) % p5_10 == 0)
            stable_10 = stable_10_p2 and stable_10_p5
            condition2_met = not stable_10
            
            if condition2_met:
                # We found the smallest n. Now we print the result and verification.
                print(f"The smallest positive integer n is: {n}")
                
                print("\n" + "="*50)
                print(f"Verification for n = {n}")
                print("="*50)

                n_minus_1 = n - 1

                print("\n--- Verifying Condition 1 (last 9 digits eventually constant) ---")
                print(f"This requires n^k(n-1) to be divisible by 10^9 = 2^9 * 5^9 for large k.")
                
                print(f"\n1a) Check for 2^9 = {p2_9}:")
                if n % 2 == 0:
                    print(f"   n = {n} is even. For k >= 9, n^k is divisible by 2^9. Condition met.")
                else:
                    q, r = divmod(n_minus_1, p2_9)
                    print(f"   n = {n} is odd. We check n-1 = {n_minus_1}.")
                    print(f"   {n_minus_1} / {p2_9} = {q} with remainder {r}.")
                    print(f"   Since the remainder is 0, n-1 is divisible by 2^9. Condition met.")
                
                print(f"\n1b) Check for 5^9 = {p5_9}:")
                if n % 5 == 0:
                    print(f"   n = {n} is a multiple of 5. For k >= 9, n^k is divisible by 5^9. Condition met.")
                else:
                    q, r = divmod(n_minus_1, p5_9)
                    print(f"   n is not a multiple of 5. We check n-1 = {n_minus_1}.")
                    print(f"   {n_minus_1} / {p5_9} = {q} with remainder {r}.")
                    print(f"   Since the remainder is 0, n-1 is divisible by 5^9. Condition met.")
                
                print("\nConclusion: Condition 1 is satisfied.")

                print("\n--- Verifying Condition 2 (last 10 digits NOT eventually constant) ---")
                print("This requires the stability condition for 10^10 to be false.")
                print("Stability for 10^10 requires BOTH of the following to be true:")
                print(f"  a) n is even OR n-1 is divisible by 2^10 = {p2_10}")
                print(f"  b) n is a multiple of 5 OR n-1 is divisible by 5^10 = {p5_10}")

                print(f"\n2a) Check for 2^10 = {p2_10}:")
                if stable_10_p2:
                     print(f"   Condition is met.")
                else:
                     q, r = divmod(n_minus_1, p2_10)
                     print(f"   n = {n} is odd and n-1 = {n_minus_1} is not divisible by {p2_10} (remainder is {r}).")
                     print(f"   Condition is NOT met.")
                
                print(f"\n2b) Check for 5^10 = {p5_10}:")
                if stable_10_p5:
                     print(f"   Condition is met.")
                else:
                     q, r = divmod(n_minus_1, p5_10)
                     print(f"   n is not a multiple of 5 and n-1 = {n_minus_1} is not divisible by {p5_10}.")
                     print(f"   Condition is NOT met.")

                print(f"\nConclusion: Since the condition for 2^10 is {'met' if stable_10_p2 else 'NOT met'} and for 5^10 is {'met' if stable_10_p5 else 'NOT met'}, their conjunction is false.")
                print("Therefore, the last 10 digits are not eventually constant. Condition 2 is satisfied.")
                
                return

    n += 1

solve()