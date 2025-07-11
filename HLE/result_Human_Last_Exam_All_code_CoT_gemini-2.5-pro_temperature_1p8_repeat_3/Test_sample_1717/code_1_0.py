def explain_max_health():
    """
    Explains the theoretical maximum health of the Eye of Cthulhu with infinite players.
    """
    base_health_expert = 3640
    per_player_increase = 0.35

    print("To find the theoretical maximum health of the Eye of Cthulhu in Expert Mode, we use the game's multiplayer health scaling formula.")
    print("-" * 40)
    
    # Explain the formula structure
    print("The formula is:")
    print("Total Health = Base Health * Health Multiplier")
    print("\nWhere the Health Multiplier is calculated as:")
    print("Health Multiplier = 1 + (Per Player Increase * (Number of Players - 1))")
    print("-" * 40)

    # Print the specific numbers for the equation
    print("For the Eye of Cthulhu in Expert Mode:")
    print(f"Base Health = {base_health_expert}")
    print(f"Per Player Increase = {per_player_increase}")
    
    print("\nSo the full equation for N players is:")
    print(f"Total Health = {base_health_expert} * (1 + {per_player_increase} * (N - 1))")
    print("-" * 40)
    
    # Explain the result with infinite players
    print("The question asks for the health given an infinite number of players (N -> ∞).")
    print("As N approaches infinity, the term '(N - 1)' also approaches infinity.")
    print("Since the health scales linearly with the number of players and there is no cap, the total health increases without bound.")
    print("\nConclusion: The theoretical maximum health is infinite.")

explain_max_health()