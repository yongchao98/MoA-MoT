import numpy as np

def find_symmetric_psne(n, u_A_func, u_B_func):
    """
    Finds all pure strategy Nash equilibria in an n-player symmetric game.

    Args:
        n (int): The number of players.
        u_A_func (function): Payoff function u(A, k) which takes k and returns the payoff.
                             k is the number of *other* players playing A.
        u_B_func (function): Payoff function u(B, k) which takes k and returns the payoff.
    """
    print(f"--- Finding PSNE for a {n}-player symmetric game ---")

    # The incentive function d(k) = u(A, k) - u(B, k)
    # It represents the incentive to play A over B when k others play A.
    # The domain for k is 0 to n-1.
    d = np.array([u_A_func(k) - u_B_func(k) for k in range(n)])
    print(f"Incentive function d(k) for k=0 to {n-1}: {d}\n")

    equilibria = []

    # A strategy profile is defined by 'm', the number of players choosing action A.
    for m in range(n + 1):
        print(f"--- Checking m = {m} ({m} players play A, {n-m} play B) ---")
        is_psne = True
        
        # Condition 1: For players choosing A (if any)
        # They should not want to switch. This requires d(m-1) >= 0.
        if m > 0:
            val = d[m-1]
            print(f"Check for A-players: Is d(m-1) >= 0?  (m-1 = {m-1})")
            print(f"Equation: d({m-1}) = {val}. Is {val} >= 0? {'Yes.' if val >= 0 else 'No.'}")
            if val < 0:
                is_psne = False
                print("Condition failed for A-players.\n")
                continue
            print("Condition met for A-players.\n")

        # Condition 2: For players choosing B (if any)
        # They should not want to switch. This requires d(m) <= 0.
        if m < n:
            val = d[m]
            print(f"Check for B-players: Is d(m) <= 0?  (m = {m})")
            print(f"Equation: d({m}) = {val}. Is {val} <= 0? {'Yes.' if val <= 0 else 'No.'}")
            if val > 0:
                is_psne = False
                print("Condition failed for B-players.\n")
                continue
            print("Condition met for B-players.\n")

        if is_psne:
            print(f"*** Profile m = {m} is a Pure Strategy Nash Equilibrium. ***\n")
            equilibria.append(m)
        
    print("--- Summary ---")
    if not equilibria:
        print("Found no Pure Strategy Nash Equilibria.")
    else:
        print(f"Found {len(equilibria)} PSNE profile(s):")
        for eq in equilibria:
            print(f"  - m = {eq} ({eq} players choose A, {n-eq} players choose B)")
    
    return equilibria


# --- Main execution: Demonstrate a game with exactly one PSNE ---
# Let's define a game that is designed to have only one equilibrium.
# We will use a simple linear incentive function.

N_PLAYERS = 5

# Let's define the payoffs such that the incentive function d(k) is strictly decreasing.
# Example: d(k) = (N_PLAYERS - 1)/2 - k
# Let's choose u(A, k) = (N_PLAYERS-1)/2 - k and u(B, k) = 0
# For n=5, u(A, k) = 2 - k
payoff_A = lambda k: (N_PLAYERS - 1) / 2.0 - k
payoff_B = lambda k: 0.0

find_symmetric_psne(N_PLAYERS, payoff_A, payoff_B)