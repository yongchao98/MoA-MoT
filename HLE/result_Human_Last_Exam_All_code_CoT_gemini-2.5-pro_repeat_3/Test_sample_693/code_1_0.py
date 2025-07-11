import sys

def calculate_revenue(bids, silent=False):
    """
    Calculates the total revenue from a divisible item auction using VCG payments.

    Args:
        bids (list): A list of integer bids.
        silent (bool): If True, suppresses print statements for internal calculations.

    Returns:
        float: The total revenue of the auction.
    """
    if not bids:
        return 0

    # Store bids with original indices to keep track of them
    indexed_bids = list(enumerate(bids))

    # --- Step 1: Determine winners with all bidders present ---
    highest_bid = 0
    if bids:
        highest_bid = max(bids)
    
    threshold = highest_bid / 2.0
    
    winners_with_all = [ib for ib in indexed_bids if ib[1] >= threshold]
    k = len(winners_with_all)

    if k == 0:
        if not silent:
            print(f"For bids {bids}, there are no winners.")
            print("Total Revenue = 0")
        return 0

    if not silent:
        print(f"--- Calculating Revenue for Bids: {bids} ---")
        winner_bids = [w[1] for w in winners_with_all]
        print(f"Highest bid is {highest_bid}. Winning threshold is {threshold}.")
        print(f"There are {k} winner(s) with bids: {winner_bids}")
        print(f"Each winner receives 1/{k} of the item.")

    # --- Step 2: Calculate VCG payment for each winner ---
    total_revenue = 0
    payment_terms = []
    
    for winner_idx, winner_bid in winners_with_all:
        # Value for others WITH this winner present
        others_value_with_winner = 0
        for other_idx, other_bid in winners_with_all:
            if other_idx != winner_idx:
                others_value_with_winner += other_bid * (1.0 / k)

        # Value for others WITHOUT this winner present
        bids_without_winner = [ib for ib in indexed_bids if ib[0] != winner_idx]
        
        others_value_without_winner = 0
        if bids_without_winner:
            other_bids_only = [b[1] for b in bids_without_winner]
            new_highest_bid = max(other_bids_only)
            new_threshold = new_highest_bid / 2.0
            
            new_winners = [ib for ib in bids_without_winner if ib[1] >= new_threshold]
            new_k = len(new_winners)
            
            if new_k > 0:
                for _, new_winner_bid in new_winners:
                    others_value_without_winner += new_winner_bid * (1.0 / new_k)

        # Payment is the difference (the "harm" caused)
        payment = others_value_without_winner - others_value_with_winner
        if not silent:
            print(f"Payment for bidder with bid {winner_bid}: {payment:.2f}")
        
        total_revenue += payment
        payment_terms.append(f"{payment:.2f}")

    if not silent:
        if len(payment_terms) > 1:
            print(f"Total Revenue = {' + '.join(payment_terms)} = {total_revenue:.2f}")
        else:
            print(f"Total Revenue = {total_revenue:.2f}")
        print("-" * 20)
        
    return total_revenue

def main():
    # Case 1: Bids are 100, 20, 5
    bids_x = [100, 20, 5]
    x = calculate_revenue(bids_x)

    # Case 2: Bids are 100, 60, 10
    bids_y = [100, 60, 10]
    y = calculate_revenue(bids_y)

    # Final output
    print(f"The final result is (x, y) = ({int(x)}, {int(y)})")

if __name__ == "__main__":
    main()