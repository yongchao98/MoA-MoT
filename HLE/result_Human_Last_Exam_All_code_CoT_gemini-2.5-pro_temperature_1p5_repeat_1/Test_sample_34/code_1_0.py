def solve_poker_question():
    """
    This function determines the worst suited jack to open from the button in a 100BB rake-free cash game.

    In poker strategy, particularly from the button, a very wide range of hands are opened to apply pressure
    on the blinds. Game Theory Optimal (GTO) charts for a 100BB rake-free game recommend opening all
    suited jack hands. These hands are JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, and J2s.

    The "worst" of these hands is the one with the lowest rank, which has the least high-card value and
    connectivity for making straights. In this list, J2s is the weakest. Therefore, it is the worst
    suited jack that you should still open from the button.
    """

    # The hand is identified based on standard GTO poker charts.
    worst_suited_jack = "J2s"

    print(f"The worst suited jack you should open from the button is: {worst_suited_jack}")

solve_poker_question()