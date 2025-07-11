def calculate_vcg_revenue(bids_with_ids, case_name):
    """
    Calculates the total revenue of the auction using VCG payment rules.
    It also prints the breakdown of the calculation.

    Args:
        bids_with_ids (list of tuples): A list of (id, bid) tuples.
        case_name (str): The name of the case ('x' or 'y') for printing.
    """
    print(f"--- Calculating revenue {case_name} for bids {[b[1] for b in bids_with_ids]} ---")

    if not bids_with_ids:
        print("No bidders, revenue is 0.")
        return 0

    # 1. Determine actual allocation
    all_bids = [b[1] for b in bids_with_ids]
    highest_bid = max(all_bids)
    threshold = 0.5 * highest_bid
    
    winners = [b for b in bids_with_ids if b[1] >= threshold]
    
    if not winners:
        print(f"Highest bid: {highest_bid}, Threshold: {threshold}. No winners.")
        print(f"The revenue {case_name} is 0.")
        return 0

    num_winners = len(winners)
    share_per_winner = 1.0 / num_winners
    
    print(f"Highest bid: {highest_bid}, Threshold: {threshold}")
    print(f"Winners (bids): {[w[1] for w in winners]}. Each gets a 1/{num_winners} share.")

    # 2. Calculate payment for each winner
    total_revenue = 0
    payment_calculations = []
    payment_values = []

    for winner_id, winner_bid in winners:
        # a. Calculate welfare of others WITH this winner present
        welfare_others_with_winner = 0
        for other_id, other_bid in winners:
            if other_id != winner_id:
                welfare_others_with_winner += other_bid * share_per_winner

        # b. Calculate welfare of others WITHOUT this winner
        bids_without_winner = [b for b in bids_with_ids if b[0] != winner_id]
        max_welfare_others_without_winner = 0
        
        if bids_without_winner:
            hypothetical_bids = [b[1] for b in bids_without_winner]
            hypo_highest_bid = max(hypothetical_bids)
            hypo_threshold = 0.5 * hypo_highest_bid
            hypo_winners = [b for b in bids_without_winner if b[1] >= hypo_threshold]

            if hypo_winners:
                hypo_num_winners = len(hypo_winners)
                hypo_share = 1.0 / hypo_num_winners
                for _, hypo_winner_bid in hypo_winners:
                    max_welfare_others_without_winner += hypo_winner_bid * hypo_share

        # c. Payment is the difference (the "harm")
        payment = max_welfare_others_without_winner - welfare_others_with_winner
        total_revenue += payment
        
        # Store for pretty printing the equation
        payment_calculations.append(f"({int(max_welfare_others_without_winner)} - {int(welfare_others_with_winner)})")
        payment_values.append(str(int(payment)))

    # 3. Print the final revenue equation
    print(f"Revenue {case_name} = {' + '.join(payment_calculations)}")
    if len(payment_values) > 1:
        print(f"         = {' + '.join(payment_values)}")
    print(f"         = {int(total_revenue)}")
    return int(total_revenue)

# Define the bids for each case. Using IDs to track bidders.
bids_x_with_ids = [(1, 100), (2, 20), (3, 5)]
bids_y_with_ids = [(1, 100), (2, 60), (3, 10)]

x = calculate_vcg_revenue(bids_x_with_ids, 'x')
print("\n" + "="*30 + "\n")
y = calculate_vcg_revenue(bids_y_with_ids, 'y')

print("\n" + "="*30 + "\n")
print(f"The final result for (x, y) is: ({x}, {y})")
