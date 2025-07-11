import sympy

def solve_score_sequences():
    """
    Calculates the number of different score sequences using generating functions.
    A match is played between team A (8 members) and team B (6 members).
    Every member of team A plays every member of team B exactly once.
    a_i = number of games won by X_i (team A)
    b_j = number of games won by Y_j (team B)
    We want to find the number of different sequences (a_1,...,a_8, b_1,...,b_6).
    """
    x = sympy.Symbol('x')

    # Parameters from the problem
    num_players_A = 8
    num_players_B = 6
    games_per_player_A = num_players_B  # an A player plays all B players
    games_per_player_B = num_players_A  # a B player plays all A players
    total_games = num_players_A * num_players_B

    # Generating function for Team A's scores
    # Each of the 8 players can win 0 to 6 games.
    # The possible scores are represented by the polynomial (1+x+...+x^6).
    # For 8 players, it's raised to the power of 8.
    poly_A = sum(x**i for i in range(games_per_player_A + 1))**num_players_A

    # Generating function for Team B's scores
    # Each of the 6 players can win 0 to 8 games.
    # The possible scores are represented by the polynomial (1+x+...+x^8).
    # For 6 players, it's raised to the power of 6.
    poly_B = sum(x**i for i in range(games_per_player_B + 1))**num_players_B

    # The generating function for the combined score sequence is the product.
    total_poly = poly_A * poly_B

    # The total number of wins must equal the total number of games.
    # We find the coefficient of x^48 in the expansion of the total polynomial.
    num_sequences = sympy.expand(total_poly).coeff(x, total_games)

    print(f"Number of players in Team A: {num_players_A}")
    print(f"Number of players in Team B: {num_players_B}")
    print(f"Games per player in Team A (max score + 1): {games_per_player_A + 1}")
    print(f"Games per player in Team B (max score + 1): {games_per_player_B + 1}")
    print(f"The total number of wins must sum to: {total_games}")
    print("\nThe number of different sequences is the coefficient of x^48 in the expansion of P(x) = (1+...+x^6)^8 * (1+...+x^8)^6.")
    print(f"Final calculated number of sequences:")
    print(num_sequences)

solve_score_sequences()