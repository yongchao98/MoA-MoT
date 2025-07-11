import math

def calculate_auction_outcome(bids_with_ids):
    """
    Helper function to determine winners and allocation value for a given set of bids.
    """
    if not bids_with_ids:
        return [], 0
    
    bid_values = [bid for _, bid in bids_with_ids]
    highest_bid = max(bid_values)
    threshold = highest_bid / 2.0
    
    winners = [bid for bid in bids_with_ids if bid[1] >= threshold]
    
    if not winners:
        return [], 0
        
    num_winners = len(winners)
    winner_bids = [bid for _, bid in winners]
    allocation_value = sum(winner_bids) / num_winners
    
    return winners, allocation_value

def calculate_revenue(bids):
    """
    Calculates the total revenue of the auction using VCG payments.
    Bids are provided as a list of numbers.
    """
    # Assign temporary unique IDs to bids for tracking
    bids_with_ids = list(enumerate(bids))

    # Determine winners in the main auction
    main_winners, _ = calculate_auction_outcome(bids_with_ids)
    
    if not main_winners:
        return 0
    
    num_main_winners = len(main_winners)
    total_revenue = 0

    # Calculate VCG payment for each winner
    for winner_id, winner_bid in main_winners:
        
        # 1. Calculate value of allocation to OTHER winners in the main auction
        other_winners_bids = [bid for w_id, bid in main_winners if w_id != winner_id]
        value_to_others_with_winner = sum(other_winners_bids) / num_main_winners
        
        # 2. Calculate value of allocation in an auction WITHOUT this winner
        bids_without_winner = [bid for bid in bids_with_ids if bid[0] != winner_id]
        _, value_of_alloc_without_winner = calculate_auction_outcome(bids_without_winner)

        # 3. Calculate VCG payment and add to revenue
        payment = value_of_alloc_without_winner - value_to_others_with_winner
        total_revenue += payment
        
    return total_revenue

def solve_auction_problem():
    """
    Solves the specific problem by calculating x and y and printing the result.
    """
    # Bids for the two scenarios
    bids_x = [100, 20, 5]
    bids_y = [100, 60, 10]

    # Calculate revenues
    x = calculate_revenue(bids_x)
    y = calculate_revenue(bids_y)

    # The final equation is (x, y) = (value_x, value_y).
    # We print the values of x and y in the specified tuple format.
    # Using math.ceil to handle potential floating point inaccuracies and round up,
    # although standard floats should be exact here. Then converting to int.
    print(f"({int(math.ceil(x - 1e-9))}, {int(math.ceil(y - 1e-9))})")

solve_auction_problem()