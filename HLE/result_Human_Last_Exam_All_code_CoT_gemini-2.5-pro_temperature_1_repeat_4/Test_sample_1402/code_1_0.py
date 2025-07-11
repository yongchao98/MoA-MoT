def solve_coin_game_probability():
    """
    Calculates the probability that Theo wins for the first time
    only after at least five games.
    """
    # From the analysis of a single game, the probability of Theo winning is 1/6.
    # The numbers in the fraction for Theo NOT winning a single game (5/6)
    p_not_win_num = 5
    p_not_win_den = 6

    # The event "first win is after at least 5 games" means "no win in the first 4 games".
    num_games_no_win = 4

    print("Step 1: Determine the probability of Theo NOT winning a single game.")
    print("The probability of Theo winning a single game is 1/6.")
    print(f"The probability of Theo NOT winning a single game is 1 - 1/6 = {p_not_win_num}/{p_not_win_den}.")
    print("\nStep 2: Calculate the probability of this happening for the first four games.")
    print("The event 'Theo wins for the first time only after at least five games' is equivalent to the event 'Theo does not win in the first four games'.")
    print("Since each game is independent, we calculate the probability as:")

    # Calculate the numerator and denominator of the final probability
    final_numerator = p_not_win_num ** num_games_no_win
    final_denominator = p_not_win_den ** num_games_no_win

    # Display the equation with all the numbers
    print(f"\n({p_not_win_num}/{p_not_win_den}) ^ {num_games_no_win} = ({p_not_win_num}^{num_games_no_win}) / ({p_not_win_den}^{num_games_no_win})")
    print(f"= {final_numerator} / {final_denominator}")

    # Calculate and display the final decimal value
    final_probability = final_numerator / final_denominator
    print(f"\nAs a decimal, the final probability is approximately: {final_probability}")

solve_coin_game_probability()