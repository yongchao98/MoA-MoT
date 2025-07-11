import sys

def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health of Terraria's
    Eye of Cthulhu in Expert Mode with an infinite number of players.
    """

    # Step 1: Define the base constants from Terraria's source code.
    base_health = 2800
    expert_mode_multiplier = 2.0
    per_player_increase = 0.35  # 35% health increase for each player after the first.

    print("To find the theoretical maximum health, we first need to establish the health formula.")
    print(f"1. The Eye of Cthulhu's base health is {base_health}.")
    print(f"2. In Expert Mode, all boss health is multiplied by {expert_mode_multiplier}.")
    print(f"3. For each player after the first, health is increased by an additional {per_player_increase*100}%.")
    print("-" * 50)

    # Step 2: Formulate the health H as a function of N players.
    # Formula: H(N) = (base_health * expert_multiplier) * (1 + (N - 1) * per_player_increase)
    # We can simplify this to a linear equation: H(N) = m*N + c

    # Calculating the slope (m) and intercept (c) of the linear equation.
    # H(N) = (5600) * (1 + 0.35N - 0.35)
    # H(N) = 5600 * (0.35N + 0.65)
    # H(N) = (5600 * 0.35) * N + (5600 * 0.65)
    slope = base_health * expert_mode_multiplier * per_player_increase
    intercept = base_health * expert_mode_multiplier * (1 - per_player_increase)

    print("The health (H) for a given number of players (N) simplifies to the linear equation: H = m*N + c")
    print("\nFinal Equation:")
    print(f"H(N) = {slope} * N + {intercept}")
    print("-" * 50)

    # Step 3: Analyze the limit as N -> infinity
    print("Because the health formula is a linear function with a positive slope, as the number")
    print("of players (N) increases, the health (H) increases without any upper limit.")
    print("\nTherefore, the theoretical maximum health for an infinite number of players is infinite.")

if __name__ == "__main__":
    calculate_eoc_health()
<<<Infinity>>>