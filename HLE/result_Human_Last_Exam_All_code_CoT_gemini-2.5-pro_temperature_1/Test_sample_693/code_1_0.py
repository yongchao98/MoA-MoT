def calculate_and_print_revenue(bids, case_name):
    """
    Calculates the revenue of a divisible item auction using VCG payments.
    """
    print(f"--- Calculating revenue '{case_name}' for bids: {bids} ---")

    if not bids:
        print("No bids provided. Revenue is 0.")
        return 0

    # Step 1: Determine winners in the original auction
    b_max = max(bids)
    threshold = b_max / 2.0
    print(f"The highest bid is {b_max}, so the winning threshold is {b_max} / 2 = {threshold}.")

    winners = []
    winner_indices = {}  # Using a dict to map original index to bid
    for i, bid in enumerate(bids):
        if bid >= threshold:
            winners.append(bid)
            winner_indices[i] = bid

    if not winners:
        print("There are no winners in this auction. Revenue is 0.")
        print(f"\nTotal revenue '{case_name}' = 0\n")
        return 0

    num_winners = len(winners)
    share = 1.0 / num_winners
    print(f"The winners are bidders with bids {winners}. There are {num_winners} winner(s), so each gets a {share:.2f} share of the item.")

    total_revenue = 0
    payments = []
    
    # Step 2: Calculate VCG payment for each winner
    for i, winner_bid in winner_indices.items():
        print(f"\nCalculating payment for the winning bidder with bid {winner_bid}:")

        # a. Calculate the welfare of OTHERS with this winner present
        welfare_others_with_winner = 0
        for j, other_winner_bid in winner_indices.items():
            if i == j:
                continue
            welfare_others_with_winner += other_winner_bid * share
        
        print(f"  Value for other winning bidders (if any) with this bidder present = {welfare_others_with_winner}")

        # b. Calculate the welfare of OTHERS in a hypothetical auction without this winner
        bids_without_winner = [b for idx, b in enumerate(bids) if idx != i]
        
        welfare_others_without_winner = 0
        if not bids_without_winner:
            print("  No other bidders remain to run a hypothetical auction.")
        else:
            print(f"  Running hypothetical auction with remaining bids: {bids_without_winner}")
            b_max_new = max(bids_without_winner)
            threshold_new = b_max_new / 2.0
            
            winners_new = [b for b in bids_without_winner if b >= threshold_new]
            
            print(f"    In this hypothetical auction, highest bid is {b_max_new}, threshold is {threshold_new}. The winners are {winners_new}.")
            
            num_winners_new = len(winners_new)
            if num_winners_new > 0:
                share_new = 1.0 / num_winners_new
                hypothetical_welfare = 0
                for winner_bid_new in winners_new:
                    hypothetical_welfare += winner_bid_new * share_new
                welfare_others_without_winner = hypothetical_welfare
        
        print(f"  Value for other bidders if this bidder were absent = {welfare_others_without_winner}")

        # c. Calculate final payment using the VCG formula
        payment = welfare_others_without_winner - welfare_others_with_winner
        print(f"  The payment is the difference: {welfare_others_without_winner} - {welfare_others_with_winner} = {payment}")
        payments.append(payment)
        total_revenue += payment
    
    payment_equation = " + ".join(map(str, payments))
    print(f"\nTotal revenue '{case_name}' is the sum of payments: {payment_equation} = {total_revenue}\n")
    return total_revenue

if __name__ == "__main__":
    # Case x
    bids_x = [100, 20, 5]
    x = calculate_and_print_revenue(bids_x, 'x')

    # Case y
    bids_y = [100, 60, 10]
    y = calculate_and_print_revenue(bids_y, 'y')

    print("-----------------------------------------")
    print(f"Final Answer: (x, y) = ({x}, {y})")
    print("-----------------------------------------")