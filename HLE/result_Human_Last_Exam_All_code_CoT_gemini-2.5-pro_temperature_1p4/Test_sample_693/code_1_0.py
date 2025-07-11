import numpy as np

def get_allocation(bids):
    """
    Calculates the allocation for a given set of bids based on the rule:
    the item is equally divided among bidders with a bid >= max_bid / 2.
    """
    if not any(b > 0 for b in bids):
        return [0.0] * len(bids)

    max_bid = float(max(bids))
    threshold = max_bid / 2.0
    
    winners = [i for i, b in enumerate(bids) if b >= threshold]
    num_winners = len(winners)
    
    allocation = [0.0] * len(bids)
    if num_winners > 0:
        share = 1.0 / num_winners
        for i in winners:
            allocation[i] = share
            
    return allocation

def calculate_payment(bidder_idx, all_bids):
    """
    Calculates the payment for a single bidder using the integral formula
    for truthful mechanisms with monotone allocation rules.
    p_i(b) = b_i * s_i(b) - integral_0^b_i s_i(z, b_-i) dz
    """
    original_bid = float(all_bids[bidder_idx])
    other_bids = [float(b) for i, b in enumerate(all_bids) if i != bidder_idx]
    
    # Current allocation for bidder i
    current_alloc = get_allocation(all_bids)[bidder_idx]
    
    # The allocation function s_i(z, b_-i) is a step function.
    # We need to find the critical points 'z' where the value of s_i(z) changes.
    critical_points = {0.0}

    # A change can occur when the bidder's status (winner/loser) changes.
    # This happens relative to the highest bid.
    
    # Case 1: Bid 'z' is NOT the highest. max_bid comes from other_bids.
    if other_bids:
        max_other_bid = max(other_bids)
        # s_i changes if z crosses the threshold max_other_bid / 2.
        critical_points.add(max_other_bid / 2.0)

    # Case 2: Bid 'z' IS the highest. max_bid = z, so threshold = z/2.
    # The number of other winners can change if z/2 crosses an other_bid b_j.
    # This happens when z = 2 * b_j.
    for b in other_bids:
        critical_points.add(2.0 * b)
        
    # We only need to integrate up to the original bid.
    sorted_points = sorted(list(critical_points))
    
    integral = 0.0
    last_point = 0.0
    
    # Integrate piecewise over the segments defined by the critical points.
    for point in sorted_points:
        if point > original_bid:
            break
        
        # Calculate s_i at a test point within the segment (last_point, point)
        if point > last_point:
            test_bid = (last_point + point) / 2.0
            temp_bids = list(other_bids)
            temp_bids.insert(bidder_idx, test_bid)
            alloc = get_allocation(temp_bids)[bidder_idx]
            integral += alloc * (point - last_point)
            
        last_point = point
        
    # Add the last segment of the integral, from the last critical point to the original_bid
    if original_bid > last_point:
        test_bid = (last_point + original_bid) / 2.0
        temp_bids = list(other_bids)
        temp_bids.insert(bidder_idx, test_bid)
        alloc = get_allocation(temp_bids)[bidder_idx]
        integral += alloc * (original_bid - last_point)

    payment = original_bid * current_alloc - integral
    return payment

def solve_and_print_revenue(name, bids):
    """
    Calculates the total revenue for a given set of bids and prints the details.
    """
    payments = [calculate_payment(i, bids) for i in range(len(bids))]
    total_revenue = sum(payments)
    
    print(f"Calculating revenue '{name}' for bids {bids}:")
    
    # Print individual payment equations
    current_allocations = get_allocation(bids)
    payment_eqs = []
    for i, p in enumerate(payments):
        # We present the numbers in the final equation as requested
        # Note: The calculation is done via the integral formula, but we show the final numbers
        payment_eqs.append(f"p_{i+1} = {p:.2f}")

    print(f"Individual payments are: {', '.join(payment_eqs)}")
    
    sum_str = " + ".join([f"{p:.2f}" for p in payments])
    print(f"Total revenue {name} = {sum_str} = {total_revenue:.2f}\n")
    return total_revenue

if __name__ == '__main__':
    # Case x: bids are 100, 20, and 5
    bids_x = [100, 20, 5]
    x = solve_and_print_revenue('x', bids_x)
    
    # Case y: bids are 100, 60, and 10
    bids_y = [100, 60, 10]
    y = solve_and_print_revenue('y', bids_y)
    
    # Final combined answer
    print(f"The final result is (x, y) = ({x:.2f}, {y:.2f})")