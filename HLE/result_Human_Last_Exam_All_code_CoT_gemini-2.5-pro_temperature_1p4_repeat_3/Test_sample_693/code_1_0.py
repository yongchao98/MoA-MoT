import sys

# Redirect stdout to a string buffer to capture print output for the final answer
# This is a trick to format the final output as requested without asking the user to copy-paste.
from io import StringIO
original_stdout = sys.stdout
sys.stdout = captured_output = StringIO()


def calculate_and_explain_revenue(bids, case_name):
    """
    Calculates the total revenue for a given set of bids in a divisible item auction.
    The function uses the VCG payment rule to ensure truthfulness.
    """
    print(f"--- Calculating revenue '{case_name}' for bids: {bids} ---")

    if not bids:
        print("No bids, revenue is 0.")
        return 0

    # 1. Determine winners and allocation with all bidders present
    b_max = max(bids)
    threshold = b_max / 2.0
    
    winners_indices = [i for i, bid in enumerate(bids) if bid >= threshold]
    
    if not winners_indices:
        print(f"The winning threshold is {threshold}. There are no winners.")
        print(f"Total revenue '{case_name}' = 0\n")
        return 0
        
    k = len(winners_indices)
    winner_bids = [bids[i] for i in winners_indices]
    print(f"Highest bid: {b_max}")
    print(f"Winning threshold (1/2 * highest bid): {threshold}")
    print(f"There are {k} winner(s) with bids: {winner_bids}. Each gets 1/{k} of the item.")

    total_revenue = 0
    payment_summary = []

    # 2. For each winner, calculate their VCG payment
    for i in winners_indices:
        winner_bid = bids[i]
        
        # a. Calculate the total value received by OTHER winners when this winner is present
        welfare_of_others_with_winner = 0
        for j in winners_indices:
            if i == j:
                continue  # Only sum over OTHER winners
            welfare_of_others_with_winner += (1.0 / k) * bids[j]

        # b. Calculate the total value others would have received if this winner were absent
        bids_without_winner = [bid for idx, bid in enumerate(bids) if idx != i]
        
        welfare_of_others_without_winner = 0
        if bids_without_winner:
            b_max_without_winner = max(bids_without_winner)
            threshold_without_winner = b_max_without_winner / 2.0
            
            # Find winners in the hypothetical scenario
            winners_indices_without_winner = [idx for idx, bid in enumerate(bids) if idx != i and bid >= threshold_without_winner]
            
            k_without_winner = len(winners_indices_without_winner)
            
            if k_without_winner > 0:
                for j_idx in winners_indices_without_winner:
                    welfare_of_others_without_winner += (1.0 / k_without_winner) * bids[j_idx]
        
        # c. The payment is the difference (the "harm" caused)
        payment_i = welfare_of_others_without_winner - welfare_of_others_with_winner
        total_revenue += payment_i
        
        print(f"\n  Calculating payment for winner with bid {winner_bid}:")
        print(f"    - Welfare of other winners (if present): {welfare_of_others_with_winner}")
        print(f"    - Welfare of other bidders (if absent): {welfare_of_others_without_winner}")
        print(f"    - Payment = {welfare_of_others_without_winner} - {welfare_of_others_with_winner} = {payment_i}")
        payment_summary.append(str(payment_i))

    print("\nFinal revenue calculation:")
    if len(payment_summary) > 1:
        print(f"Total Revenue '{case_name}' = {' + '.join(payment_summary)} = {total_revenue}")
    else:
        print(f"Total Revenue '{case_name}' = {total_revenue}")
    print("-" * (len(case_name) + 30))
    print() # Add a newline for separation
    return total_revenue

# Bids for case x
bids_x = [100, 20, 5]
x = calculate_and_explain_revenue(bids_x, "x")

# Bids for case y
bids_y = [100, 60, 10]
y = calculate_and_explain_revenue(bids_y, "y")

# Restore original stdout
sys.stdout = original_stdout
# Print the captured output to the console
captured_text = captured_output.getvalue()
print(captured_text)

# Final answer
final_answer_tuple = (x, y)
print(f"The final result for (x, y) is: {final_answer_tuple}")
print(f"<<<{final_answer_tuple}>>>")
