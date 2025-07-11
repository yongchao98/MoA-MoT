def calculate_and_print_revenue(bids, variable_name):
    """
    Calculates and prints the revenue for the divisible item auction.
    
    Args:
        bids (list): A list of numerical bids.
        variable_name (str): The name of the revenue variable ('x' or 'y').
    """
    print(f"--- Calculating revenue '{variable_name}' for bids: {bids} ---")
    
    if not bids:
        print("No bids provided. Revenue is 0.")
        return 0

    highest_bid = max(bids)
    winning_threshold = highest_bid / 2

    print(f"The highest bid is {highest_bid}.")
    print(f"The winning threshold is {highest_bid} / 2 = {winning_threshold}.")

    winners_info = []
    for i, bid in enumerate(bids):
        if bid >= winning_threshold:
            other_bids = bids[:i] + bids[i+1:]
            highest_other_bid = max(other_bids) if other_bids else 0
            payment = highest_other_bid / 2
            winners_info.append({
                'bid': bid,
                'payment': payment,
                'highest_other_bid': highest_other_bid
            })

    if not winners_info:
        print("There are no winners.")
        total_revenue = 0
        print(f"Total revenue {variable_name} = 0")
        return total_revenue

    print("\nIndividual payment calculations for winners:")
    payments = []
    for winner in winners_info:
        payments.append(winner['payment'])
        print(f"  - Payment for bid {winner['bid']}: max of other bids is {winner['highest_other_bid']}, so payment is {winner['highest_other_bid']} / 2 = {winner['payment']}.")

    total_revenue = sum(payments)
    
    # Constructing the final equation string
    equation_str = " + ".join(map(lambda p: str(int(p)) if p.is_integer() else str(p), payments))
    if len(payments) > 1:
        final_revenue_str = str(int(total_revenue)) if total_revenue.is_integer() else str(total_revenue)
        print(f"\nTotal revenue {variable_name} = {equation_str} = {final_revenue_str}")
    else:
        final_revenue_str = str(int(total_revenue)) if total_revenue.is_integer() else str(total_revenue)
        print(f"\nTotal revenue {variable_name} = {final_revenue_str}")

    return total_revenue

# Calculate x for the first set of bids
bids_x = [100, 20, 5]
x = calculate_and_print_revenue(bids_x, 'x')

print("\n" + "="*50 + "\n")

# Calculate y for the second set of bids
bids_y = [100, 60, 10]
y = calculate_and_print_revenue(bids_y, 'y')

print("\n" + "="*50)
print(f"The final result is (x, y) = ({int(x)}, {int(y)})")
<<< (10, 80) >>>