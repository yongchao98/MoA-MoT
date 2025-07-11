import collections

def get_allocations(bids_with_indices):
    """
    Determines the allocation for each bidder.
    Returns a dictionary mapping bidder index to their share of the item.
    """
    if not bids_with_indices:
        return {}

    bids = [b for _, b in bids_with_indices]
    highest_bid = max(bids)
    threshold = highest_bid / 2.0
    
    winner_indices = [i for i, b in bids_with_indices if b >= threshold]
    
    if not winner_indices:
        return {}
        
    num_winners = len(winner_indices)
    share = 1.0 / num_winners
    
    allocations = {idx: share for idx in winner_indices}
    return allocations

def calculate_revenue(bids, case_name):
    """
    Calculates the total revenue for a given set of bids using VCG payments.
    """
    print(f"--- Calculating revenue '{case_name}' for bids: {bids} ---")
    
    bids_with_indices = list(enumerate(bids))
    
    # 1. Determine allocations with everyone present
    allocations_with_all = get_allocations(bids_with_indices)
    
    total_revenue = 0
    payments = []
    
    for i, bid_i in bids_with_indices:
        # Losers pay 0
        if i not in allocations_with_all:
            payments.append(0)
            print(f"Bidder with bid {bid_i} is a loser. Payment is 0.")
            continue
        
        # This is a winner, so we calculate their VCG payment
        
        # a) Calculate total value for others if bidder 'i' is absent
        bids_without_i = [(idx, b) for idx, b in bids_with_indices if idx != i]
        # We need to re-index for the helper function to work correctly
        reindexed_bids_without_i = list(enumerate([b for _, b in bids_without_i]))
        original_indices_map = {new_idx: orig_idx for new_idx, (orig_idx, _) in enumerate(bids_without_i)}
        
        allocs_without_i = get_allocations(reindexed_bids_without_i)
        
        sw_others_without_i = 0
        for new_idx, alloc_share in allocs_without_i.items():
            original_bidder_index = original_indices_map[new_idx]
            original_bid_value = bids[original_bidder_index]
            sw_others_without_i += original_bid_value * alloc_share
            
        # b) Calculate total value for others if bidder 'i' is present
        sw_others_with_i = 0
        for j, bid_j in bids_with_indices:
            if i == j:
                continue
            sw_others_with_i += bid_j * allocations_with_all.get(j, 0)
            
        # c) Payment is the difference (the "harm" caused)
        payment_i = sw_others_without_i - sw_others_with_i
        payments.append(payment_i)
        total_revenue += payment_i
        
        print(f"Bidder with bid {bid_i} is a winner.")
        print(f"  - Social welfare of others if this bidder is absent: {sw_others_without_i:.2f}")
        print(f"  - Social welfare of others if this bidder is present: {sw_others_with_i:.2f}")
        print(f"  - Payment = {sw_others_without_i:.2f} - {sw_others_with_i:.2f} = {payment_i:.2f}")

    print("\nFinal revenue equation:")
    payment_str = " + ".join([f"{p:.2f}" for p in payments])
    print(f"Total revenue for '{case_name}' = {payment_str} = {total_revenue:.2f}")
    print("-" * (len(case_name) + 33))
    
    return total_revenue

if __name__ == "__main__":
    bids_x = [100, 20, 5]
    x = calculate_revenue(bids_x, 'x')
    
    print("\n")
    
    bids_y = [100, 60, 10]
    y = calculate_revenue(bids_y, 'y')
    
    print("\nFinal Result:")
    # The output format (x, y) as requested
    final_result = (x, y)
    print(f"(x, y) = {final_result}")
