import math

def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health for the Eye of Cthulhu
    in Expert Mode with an infinite number of players.
    """
    # Step 1: Define the base parameters for the Eye of Cthulhu in Expert Mode.
    base_health = 3640
    player_scaling_factor = 0.35

    # Step 2: Explain the health scaling formula.
    # The formula is: Total Health = BaseHealth * (1 + ScalingFactor * (NumberOfPlayers - 1))
    # Let N be the number of players. The formula becomes: H(N) = 3640 * (1 + 0.35 * (N - 1))
    # This simplifies to H(N) = 3640 * (1 + 0.35*N - 0.35) = 3640 * (0.65 + 0.35*N).

    # Step 3: Analyze the formula for an infinite number of players.
    # As the number of players (N) approaches infinity, the term (0.65 + 0.35*N) also
    # approaches infinity. This means the total health increases without any upper bound.
    # Therefore, the theoretical maximum health is infinite.

    # Step 4: Print the final equation with the symbolic representation.
    num_players_symbol = "N (approaching infinity)"
    final_result = "Infinity"

    print("Terraria Expert Mode Health Calculation for Eye of Cthulhu")
    print("-" * 60)
    print(f"Base Health: {base_health}")
    print(f"Per-Player Health Scaling Factor: {player_scaling_factor}")
    print("\nFormula: Base Health * (1 + Scaling Factor * (Number of Players - 1))")
    print("\nAs the Number of Players approaches infinity, the equation becomes:")

    # We use strings to represent the symbolic nature of the calculation.
    equation_str = f"{base_health} * (1 + {player_scaling_factor} * ({num_players_symbol} - 1)) = {final_result}"
    print(equation_str)

if __name__ == "__main__":
    calculate_eoc_health()