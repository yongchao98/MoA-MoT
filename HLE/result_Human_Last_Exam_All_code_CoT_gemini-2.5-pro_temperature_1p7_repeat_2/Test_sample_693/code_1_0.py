def solve_auction():
    """
    Calculates auction revenue based on a truthful payment rule.

    The allocation rule is that an item is equally divided among all bidders
    whose bid is at least half the highest bid. For the auction to be truthful,
    we must define a payment rule. The "highest losing bid" rule is a plausible
    candidate that ensures individual rationality (non-negative utility for truthful bidders).
    Under this rule, every winner pays a price equal to the highest bid among the losers.
    """

    def calculate_revenue(bids, scenario_name):
        """
        Calculates the total revenue for a given set of bids.
        """
        print(f"--- Calculating Revenue for Scenario {scenario_name} ---")
        print(f"The bids are: {bids}")

        if not bids:
            print("No bids, so revenue is 0.")
            return 0

        highest_bid = max(bids)
        threshold = 0.5 * highest_bid

        winners = []
        losers = []
        for bid in bids:
            if bid >= threshold:
                winners.append(bid)
            else:
                losers.append(bid)

        print(f"The highest bid is {highest_bid}.")
        print(f"The winning threshold is 1/2 * {highest_bid} = {threshold}.")
        print(f"Winning bids: {winners}")
        print(f"Losing bids: {losers}")

        if not winners:
            print("No winners, so revenue is 0.")
            return 0
        
        payment_per_winner = 0
        if losers:
            payment_per_winner = max(losers)
            print(f"The payment for each winner is the highest losing bid, which is {payment_per_winner}.")
        else:
            print("No losers, so payment is 0.")
            
        num_winners = len(winners)
        total_revenue = num_winners * payment_per_winner

        print(f"There are {num_winners} winner(s).")
        print(f"The total revenue {scenario_name} is {num_winners} * {payment_per_winner} = {total_revenue}.\n")
        return total_revenue

    # Scenario for x
    bids_x = [100, 20, 5]
    x = calculate_revenue(bids_x, 'x')

    # Scenario for y
    bids_y = [100, 60, 10]
    y = calculate_revenue(bids_y, 'y')

    print(f"Final Answer: (x, y) = ({x}, {y})")


solve_auction()