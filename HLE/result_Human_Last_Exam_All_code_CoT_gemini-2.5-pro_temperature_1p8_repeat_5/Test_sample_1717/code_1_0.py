import math

def calculate_eoc_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """

    # Step 1: Define base parameters for Expert Mode Eye of Cthulhu
    base_health_expert = 3640
    
    # Step 2: Define the multiplayer health scaling formula.
    # In Expert Mode, health scales by an additional 35% for each player after the first.
    # Formula: TotalHealth = BaseHealth * (1 + 0.35 * (NumberOfPlayers - 1))
    # Let N be the NumberOfPlayers.
    # H(N) = 3640 * (1 + 0.35 * (N - 1))
    
    # Step 3: Simplify the formula into a linear equation: H(N) = m*N + b
    # H(N) = 3640 * (1 + 0.35*N - 0.35)
    # H(N) = 3640 * (0.65 + 0.35*N)
    # H(N) = (3640 * 0.35) * N + (3640 * 0.65)
    
    # Calculate the coefficients 'm' (per-player health) and 'b' (the constant base)
    m = base_health_expert * 0.35
    b = base_health_expert * (1 - 0.35)
    
    print("Terraria Expert Mode Health Scaling Formula:")
    print("Total Health = Base Health * (1 + 0.35 * (Number of Players - 1))")
    print(f"For the Eye of Cthulhu, the Base Health is {base_health_expert}.")
    print("\nSimplifying this for the Eye of Cthulhu gives a linear equation:")
    print("Total Health = (m * Number of Players) + b")
    print("\nWhere the numbers in the final equation are:")
    print(f"Health increase per player (m): {m}")
    print(f"Adjusted base health (b): {b}")
    
    # Step 4: Analyze the result for an infinite number of players.
    print("\nThe resulting health formula is: Health =", m, "* (Number of Players) +", b)
    print("\nAs you can see, the health is a linear function of the number of players.")
    print("This means that as the number of players (N) increases, the health increases without any upper limit.")
    print("\nTherefore, as the number of players approaches infinity, the theoretical health also approaches infinity.")

calculate_eoc_health()