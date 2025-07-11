def solve_auction_revenue():
    """
    Calculates the auction revenue for two scenarios based on a specific divisible item auction rule.
    """
    # Case 1: Bids for calculating x
    bids_x = [100, 20, 5]
    highest_bid_x = max(bids_x)
    
    # In a truthful auction with the given rules, the total revenue
    # is 1/2 times the highest bid.
    revenue_x = 0.5 * highest_bid_x
    
    print("Calculating revenue x:")
    print(f"The bids are: {bids_x}")
    print(f"The highest bid is: {highest_bid_x}")
    print(f"The revenue x is calculated as: 0.5 * {highest_bid_x} = {revenue_x}")
    print("-" * 30)
    
    # Case 2: Bids for calculating y
    bids_y = [100, 60, 10]
    highest_bid_y = max(bids_y)
    
    # The revenue rule remains the same.
    revenue_y = 0.5 * highest_bid_y
    
    print("Calculating revenue y:")
    print(f"The bids are: {bids_y}")
    print(f"The highest bid is: {highest_bid_y}")
    print(f"The revenue y is calculated as: 0.5 * {highest_bid_y} = {revenue_y}")
    print("-" * 30)
    
    # Final result
    print(f"The final result for (x, y) is: ({int(revenue_x)}, {int(revenue_y)})")

solve_auction_revenue()