def solve_bike_race_problem():
    """
    This function prints the derived equations for the achieved growth rate W
    and the decrease in growth rate ΔW.
    The derivation is explained in the text preceding this code block.
    """

    # The achieved doubling rate W is given by the sum of p_i * log(Return_i)
    # where p_i are the true probabilities and Return_i are the wealth multipliers
    # from betting based on the incorrect probabilities.
    # p = [1/2, 1/4, 1/8, 1/8]
    # Bet: 1/4 on bike 2.
    # Return_1 = 3/4, Return_2 = 3/2, Return_3 = 3/4, Return_4 = 3/4
    w_equation = "W = (1/2) * log(3/4) + (1/4) * log(3/2) + (1/8) * log(3/4) + (1/8) * log(3/4)"

    # The decrease in the doubling rate ΔW is W* - W.
    # W* = log(2) - (1/2)log(3)
    # W  = log(3) - (7/4)log(2)
    # ΔW = W* - W = (11/4)log(2) - (3/2)log(3)
    delta_w_equation = "ΔW = (11/4) * log(2) - (3/2) * log(3)"

    print("The doubling rate you will achieve is:")
    print(w_equation)
    print("\nThe decrease in your doubling rate is:")
    print(delta_w_equation)

solve_bike_race_problem()