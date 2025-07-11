import sympy

def calculate_eoc_max_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Expert Mode
    with an infinite number of players.
    """

    # According to the Terraria Wiki, which documents information from the game's source code:
    # 1. The base health for the Eye of Cthulhu in Expert Mode is 3640.
    # 2. In Expert Mode, boss health increases by 35% (0.35) for each player beyond the first.

    base_health = 3640
    life_bonus_per_player = 0.35

    # We define the number of players as a mathematical symbol 'N'.
    N = sympy.Symbol('N')

    # The formula for health scaling in multiplayer Expert Mode is:
    # Health = BaseHealth * (1 + LifeBonusPerPlayer * (NumberOfPlayers - 1))
    health_equation = base_health * (1 + life_bonus_per_player * (N - 1))

    # To find the theoretical maximum with infinite players, we calculate the limit
    # of the health equation as N approaches infinity (oo).
    theoretical_max_health = sympy.limit(health_equation, N, sympy.oo)

    # Let's simplify the equation to see the linear relationship.
    # H = 3640 * (1 + 0.35*N - 0.35)
    # H = 3640 * (0.65 + 0.35*N)
    # H = (3640 * 0.65) + (3640 * 0.35) * N
    # H = 2366 + 1274 * N
    constant_part = base_health * (1 - life_bonus_per_player)
    scaling_part = base_health * life_bonus_per_player

    print("Terraria Boss Health Calculation: Eye of Cthulhu (Expert Mode)")
    print("-" * 60)
    print(f"Base Health: {base_health}")
    print(f"Life Bonus Per Additional Player: {life_bonus_per_player}")
    print("\nThe general formula is: H = BaseHealth * (1 + Bonus * (N-1))")
    print("\nSimplified linear equation:")
    print(f"H = {int(constant_part)} + {int(scaling_part)} * N")
    print("\nAs the number of players (N) approaches infinity, the total health also approaches infinity.")
    print(f"\nThe theoretical maximum health is: {theoretical_max_health}")

if __name__ == '__main__':
    calculate_eoc_max_health()
