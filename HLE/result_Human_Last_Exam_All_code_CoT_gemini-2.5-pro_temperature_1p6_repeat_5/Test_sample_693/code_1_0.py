def calculate_auction_revenue(bids, label):
    """
    Calculates the total revenue from a divisible item auction.

    The auction rules are:
    1. Allocation: Item is equally divided among bidders with bid >= 0.5 * highest_bid.
    2. Payment: The auction is truthful. We assume this implies winners pay their
       critical price, which is the lowest bid they needed to win.
       This critical price for bidder i is 0.5 * max(bids of others).
    """
    print(f"Calculating revenue '{label}' for bids: {bids}")

    if not bids:
        print("No bids, revenue is 0.")
        return 0

    highest_bid = max(bids)
    win_threshold = 0.5 * highest_bid
    print(f"The highest bid is {highest_bid}.")
    print(f"The winning threshold is 0.5 * {highest_bid} = {win_threshold}")

    winners = []
    winner_payments = []
    
    # Identify winners and their payments
    for i, bid in enumerate(bids):
        if bid >= win_threshold:
            other_bids = bids[:i] + bids[i+1:]
            if not other_bids:
                # If only one bidder, they pay 0.
                highest_other_bid = 0
            else:
                highest_other_bid = max(other_bids)
            
            payment = 0.5 * highest_other_bid
            winners.append(bid)
            winner_payments.append(payment)
            print(f"Bidder with bid {bid} is a WINNER.")
            print(f"  - Payment calculation: 0.5 * max(other bids) = 0.5 * {highest_other_bid} = {payment}")
        else:
            print(f"Bidder with bid {bid} is a LOSER.")

    # Calculate total revenue
    total_revenue = sum(winner_payments)

    if winner_payments:
        payment_summary = " + ".join(map(str, winner_payments))
        print(f"Total revenue '{label}' is the sum of payments: {payment_summary} = {total_revenue}")
    else:
        print(f"No winners, so total revenue '{label}' is 0.")

    print("-" * 20)
    return total_revenue

if __name__ == "__main__":
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]

    x = calculate_auction_revenue(bids_x, 'x')
    y = calculate_auction_revenue(bids_y, 'y')

    print(f"Final result (x, y) = ({x}, {y})")
    # The final answer format is specified by the system.
    # The output format below is to match the expected answer format.
    final_answer_tuple = (x, y)
    print(f"<<<{final_answer_tuple}>>>")
