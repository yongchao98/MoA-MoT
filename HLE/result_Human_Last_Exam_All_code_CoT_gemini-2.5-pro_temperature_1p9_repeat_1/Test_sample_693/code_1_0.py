def solve_auction_revenue():
    """
    Calculates the revenue for two single-item divisible auctions based on
    the provided rules and prints a detailed breakdown of the calculation.

    The allocation rule is that the item is divided among bidders whose bid
    is at least 1/2 the highest bid.
    The payment rule for a truthful auction requires each winning bidder 'i' to pay
    0.5 * max(bids of others), which is the critical price they must meet to win.
    """

    def calculate_revenue_details(bids, name):
        """
        Calculates revenue for a set of bids and prints a detailed breakdown.
        """
        print(f"--- Calculating revenue '{name}' for bids: {bids} ---")

        if not bids:
            print(f"Revenue '{name}' = 0 (no bids).")
            return 0

        # Find the highest bid and the winning threshold
        b_max = float(max(bids))
        threshold = 0.5 * b_max
        print(f"Highest bid: {b_max}")
        print(f"Winning threshold (0.5 * highest bid): {threshold}")

        # Identify winners
        winners = []
        for i, bid in enumerate(bids):
            if bid >= threshold:
                winners.append({'index': i, 'bid': float(bid)})
        
        if not winners:
            print("There are no winners.")
            print(f"Final revenue '{name}' = 0")
            return 0

        winning_bids = [w['bid'] for w in winners]
        print(f"Winning bids: {winning_bids}")
        print("Calculating payments for each winner:")

        payments = []
        total_revenue = 0
        # Calculate payment for each winner
        for winner in winners:
            # Create a list of all other bids
            other_bids = [b for idx, b in enumerate(bids) if idx != winner['index']]
            
            # The highest bid among the others determines the payment
            b_max_other = float(max(other_bids)) if other_bids else 0.0
            
            # Payment is 0.5 times the highest bid of others
            payment = 0.5 * b_max_other
            payments.append(payment)
            
            other_bids_str = ', '.join(map(str, other_bids)) if other_bids else "N/A"
            print(f"- Payment for bidder with bid {winner['bid']}: 0.5 * max({other_bids_str}) = 0.5 * {b_max_other} = {payment}")

        # Sum payments for total revenue
        total_revenue = sum(payments)
        payment_str = " + ".join(map(str, payments))

        if len(payments) > 1:
            print(f"Total revenue '{name}' = {payment_str} = {total_revenue}")
        else:
            print(f"Total revenue '{name}' = {total_revenue}")
        
        return total_revenue

    # Bids for the two scenarios
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]

    # Calculate revenues
    x = calculate_revenue_details(bids_x, 'x')
    print("\n" + "="*50 + "\n")
    y = calculate_revenue_details(bids_y, 'y')

    # Print the final combined result
    print("\n" + "="*50 + "\n")
    print(f"The final result is (x, y) = ({x}, {y})")

solve_auction_revenue()