import sys
import io

# Redirect stdout to capture prints
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def get_allocation(bids_with_indices):
    """
    Calculates the allocation of the item based on the bids.
    Takes a list of (original_index, bid) tuples.
    Returns a dictionary mapping original_index to share, and a list of winner original_indices.
    """
    if not bids_with_indices:
        return {}, []

    # Find the highest bid
    highest_bid = 0
    if bids_with_indices:
        highest_bid = max(bids_with_indices, key=lambda item: item[1])[1]

    if highest_bid == 0:
        allocation = {original_index: 0 for original_index, bid in bids_with_indices}
        return allocation, []

    # Determine winners
    threshold = highest_bid / 2.0
    winners = []
    for original_index, bid in bids_with_indices:
        if bid >= threshold:
            winners.append(original_index)

    # Calculate allocation
    num_winners = len(winners)
    allocation = {original_index: 0 for original_index, bid in bids_with_indices}
    if num_winners > 0:
        share = 1.0 / num_winners
        for winner_index in winners:
            allocation[winner_index] = share
            
    return allocation, winners

def calculate_and_explain_revenue(scenario_name, bids):
    """
    Calculates the VCG revenue for a given set of bids and explains the process.
    """
    print(f"### Calculation for Revenue '{scenario_name}' ###")
    print(f"Bids: {bids}")

    # The bids are treated as the bidders' true valuations
    valuations = bids
    
    # Create a list of (original_index, bid) tuples to keep track of bidders
    bids_with_indices = list(enumerate(valuations))
    
    # 1. Calculate the outcome with all bidders present
    base_allocation, base_winners = get_allocation(bids_with_indices)
    
    if not base_winners:
        print("There are no winners.")
        print(f"The revenue {scenario_name} is 0.\n")
        return 0, []

    print(f"Highest bid is {max(valuations)}. Winning threshold is {max(valuations)/2.0}.")
    winner_bids = [bids[i] for i in sorted(base_winners)]
    print(f"There are {len(base_winners)} winner(s) with bids: {winner_bids}")
    
    total_revenue = 0
    payments = []
    payment_dict = {}
    
    # 2. Calculate payment for each winner based on VCG
    for winner_index in sorted(base_winners):
        winner_valuation = valuations[winner_index]
        
        # Calculate welfare of OTHER bidders WITH the current winner present
        welfare_others_with_winner = 0
        for i, val in enumerate(valuations):
            if i != winner_index:
                welfare_others_with_winner += val * base_allocation.get(i, 0)

        # Calculate welfare of OTHER bidders WITHOUT the current winner
        bids_without_winner = [(idx, bid) for idx, bid in bids_with_indices if idx != winner_index]
        counterfactual_alloc, _ = get_allocation(bids_without_winner)
        
        welfare_others_without_winner = 0
        for idx, share in counterfactual_alloc.items():
             welfare_others_without_winner += valuations[idx] * share

        # VCG Payment is the "harm" caused to others
        payment = welfare_others_without_winner - welfare_others_with_winner
        payments.append(payment)
        payment_dict[winner_index] = payment
        total_revenue += payment
    
    print("\nCalculating payments for winners:")
    for winner_index in sorted(base_winners):
        print(f"  - Payment for bidder with bid {bids[winner_index]}: {payment_dict[winner_index]}")

    # 3. Print the final revenue equation
    payment_strings = [str(round(p, 2)) for p in payments]
    if not payment_strings:
        equation = f"{scenario_name} = 0"
    else:
        equation = f"{scenario_name} = {' + '.join(payment_strings)}"
    
    if len(payment_strings) > 1:
        equation += f" = {round(total_revenue, 2)}"
        
    print(f"\nTotal Revenue Equation: {equation}")
    print("-" * 20)
    return total_revenue, payments

# --- Main execution ---
bids_x = [100, 20, 5]
bids_y = [100, 60, 10]

x, _ = calculate_and_explain_revenue('x', bids_x)
y, _ = calculate_and_explain_revenue('y', bids_y)

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

final_answer_tuple = (round(x), round(y))
print(f"The final result for (x, y) is: {final_answer_tuple}")

# Final Answer Format
print(f"<<<{final_answer_tuple}>>>")
