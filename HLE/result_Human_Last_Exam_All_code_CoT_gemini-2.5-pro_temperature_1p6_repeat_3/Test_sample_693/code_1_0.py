def calculate_auction_revenue(bids, name):
    """
    Calculates the revenue of a divisible item auction based on the specified rules.

    Args:
        bids (list of int/float): The list of bids from all bidders.
        name (str): The name of the revenue variable (e.g., 'x' or 'y').

    Returns:
        float: The total revenue from the auction.
    """
    if not bids:
        print(f"Revenue {name} = 0 (No bids)")
        return 0

    # Find the highest bid
    highest_bid = max(bids)

    # Determine the winning threshold price
    threshold_price = 0.5 * highest_bid

    # Count the number of bidders whose bid is at least the threshold
    num_winners = sum(1 for bid in bids if bid >= threshold_price)

    # Calculate the total revenue
    revenue = num_winners * threshold_price
    
    # Print the equation as requested
    print(f"For bids {bids}:")
    print(f"The winning threshold is 1/2 * {highest_bid} = {threshold_price}")
    print(f"The number of winners is {num_winners}")
    print(f"Revenue {name} = {num_winners} * {threshold_price} = {revenue}")
    
    return revenue

# Bids for the first scenario
bids_x = [100, 20, 5]
x = calculate_auction_revenue(bids_x, 'x')

print("-" * 20)

# Bids for the second scenario
bids_y = [100, 60, 10]
y = calculate_auction_revenue(bids_y, 'y')

print("-" * 20)

# Print the final result tuple (x, y)
# Using int() to display clean integers if they have no fractional part.
print(f"The final result (x, y) is: ({int(x)}, {int(y)})")
