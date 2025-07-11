def calculate_auction_revenue(bids, case_name):
    """
    Calculates the total revenue from the described divisible item auction.

    Args:
        bids (list): A list of numbers representing the bids.
        case_name (str): The name of the variable for the revenue (e.g., 'x' or 'y').

    Returns:
        tuple: A tuple containing the total revenue and a list of individual payments.
    """
    if not bids:
        return 0, []

    # 1. Find the highest bid and the winning threshold
    b_max = max(bids)
    threshold = b_max / 2

    # 2. Identify the winners
    winners = []
    for bid in bids:
        if bid >= threshold:
            winners.append(bid)

    # 3. For each winner, calculate their payment and sum for total revenue
    total_revenue = 0
    payments = []
    if not winners:
        print(f"For bids {bids}:")
        print(f"Highest bid is {b_max}, threshold is {threshold}. There are no winners.")
        print(f"{case_name} = 0")
        return 0, []

    for winner_bid in sorted(winners, reverse=True):
        # Create a list of the other bids
        other_bids = list(bids)
        other_bids.remove(winner_bid)
        
        # The payment is half the highest bid among the others
        b_max_other = max(other_bids) if other_bids else 0
        payment = b_max_other / 2
        payments.append(payment)
        total_revenue += payment

    # 4. Print the detailed results
    print(f"For bids {bids}:")
    print(f"Highest bid is {b_max}, winning threshold is {b_max / 2:.2f}. Winners are bidders with bids: {sorted(winners, reverse=True)}.")
    
    # Create the equation string
    payment_str = " + ".join([f"{p:.0f}" for p in sorted(payments, reverse=True)])
    print(f"The payments for the winners are: {sorted([p for p in payments], reverse=True)}")
    print(f"{case_name} = {payment_str} = {total_revenue:.0f}")
    
    return total_revenue, payments

# --- Main execution ---

# Case 1: Calculate x
bids_x = [100, 20, 5]
x, _ = calculate_auction_revenue(bids_x, 'x')

print("\n" + "-"*40 + "\n")

# Case 2: Calculate y
bids_y = [100, 60, 10]
y, _ = calculate_auction_revenue(bids_y, 'y')

# Final Answer
print("\n" + "-"*40 + "\n")
print(f"The final result is (x, y) = ({x:.0f}, {y:.0f})")

final_answer = (int(x), int(y))
# The final answer is formatted as requested below.
# <<< (10, 80) >>>