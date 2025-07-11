def solve_auction(bids, name):
    """
    Calculates the auction revenue based on the described mechanism.
    A truthful payment rule for this allocation results in a total
    revenue equal to half of the highest bid.
    """
    if not bids:
        highest_bid = 0
    else:
        highest_bid = max(bids)

    revenue = highest_bid / 2

    # Output the equation for the revenue calculation as requested.
    # Note: Using int() for clean display as results are whole numbers.
    print(f"Revenue {name} = {highest_bid} / 2 = {int(revenue)}")
    return int(revenue)

# Calculate revenue x for the first set of bids
print("Calculating revenue x for bids [100, 20, 5]:")
x = solve_auction([100, 20, 5], 'x')

# Calculate revenue y for the second set of bids
print("\nCalculating revenue y for bids [100, 60, 10]:")
y = solve_auction([100, 60, 10], 'y')

# Print the final result tuple
print(f"\nThe pair (x, y) is: ({x}, {y})")
