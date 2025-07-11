def find_symmetric_psne(n, payoffs):
    """
    Finds the Pure Strategy Nash Equilibria (PSNE) in an n-player symmetric game
    by analyzing the incentive to play action 'A' over 'B'.

    Args:
        n (int): The number of players.
        payoffs (dict): A dictionary with keys 'A' and 'B'.
                        payoffs['A'][k] is the payoff for action 'A' when k others play 'A'.
                        payoffs['B'][k] is the payoff for action 'B' when k others play 'A'.
                        Both should be lists of length n.
    """
    if len(payoffs['A']) != n or len(payoffs['B']) != n:
        raise ValueError(f"Payoff vectors must have length n={n}")

    # d(k) = P(A, k) - P(B, k) represents the incentive for a player to
    # choose action 'A' over 'B', given k other players chose 'A'.
    d = [payoffs['A'][k] - payoffs['B'][k] for k in range(n)]
    
    print(f"Analyzing an n={n} player game.")
    print(f"Payoff function for A, P(A,k): {payoffs['A']}")
    print(f"Payoff function for B, P(B,k): {payoffs['B']}")
    print(f"Incentive function d(k) = P(A,k) - P(B,k): {[f'{x:.2f}' for x in d]}\n")
    
    psne_count = 0
    
    # Check for PSNE where m=0 players choose 'A' (all 'B')
    # Condition: A player choosing 'B' must not want to switch.
    # The incentive to switch to 'A' is d(0). We need d(0) <= 0.
    if d[0] <= 0:
        psne_count += 1
        print(f"Found PSNE: m = 0 players choose A (all B).")
        # Final equation output per instruction
        print(f"--> Condition check: d(0) = {d[0]:.2f}, which is <= 0.\n")
        
    # Check for PSNE for m in {1, ..., n-1}
    for m in range(1, n):
        # Cond 1: A player choosing 'A' must not want to switch to 'B' (d(m-1) >= 0)
        # Cond 2: A player choosing 'B' must not want to switch to 'A' (d(m) <= 0)
        cond1 = d[m-1] >= 0
        cond2 = d[m] <= 0
        if cond1 and cond2:
            psne_count += 1
            print(f"Found PSNE: m = {m} players choose A.")
            # Final equation output per instruction
            print(f"--> Condition check: d({m-1}) = {d[m-1]:.2f} (>= 0) and d({m}) = {d[m]:.2f} (<= 0).\n")

    # Check for PSNE where m=n players choose 'A' (all 'A')
    # Condition: A player choosing 'A' must not want to switch.
    # The incentive to prefer 'A' is d(n-1). We need d(n-1) >= 0.
    if d[n-1] >= 0:
        psne_count += 1
        print(f"Found PSNE: m = n = {n} players choose A (all A).")
        # Final equation output per instruction
        print(f"--> Condition check: d({n-1}) = {d[n-1]:.2f}, which is >= 0.\n")

    print(f"Total number of symmetric equilibrium types found: {psne_count}")


# --- Demonstration ---
# We define a game for n=5 players that we expect to have exactly one PSNE.
# In this game, action 'A' is always preferable to 'B'.
n_players = 5
game_payoffs = {
    'A': [1.0, 1.0, 1.0, 1.0, 1.0],  # Payoff for 'A' is always 1
    'B': [0.0, 0.0, 0.0, 0.0, 0.0]   # Payoff for 'B' is always 0
}

find_symmetric_psne(n_players, game_payoffs)
