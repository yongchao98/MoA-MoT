import math

def calculate_returns(trade_price, subsequent_mid_price):
    """
    Calculates and demonstrates the 15s return for a buyer and a seller.

    The core finding is that the buyer's return is higher than the seller's,
    which implies the subsequent price p' is higher than the trade price p.
    This code demonstrates this calculation with an example.
    """
    p = trade_price
    p_prime = subsequent_mid_price

    print(f"Scenario Demonstration:")
    print(f"Initial Trade Price (p): {p}")
    print(f"Mid-Price 15s later (p'): {p_prime}\n")

    # Calculate the base return value (R = p'/p - 1)
    base_return = (p_prime / p) - 1

    # Buyer's return is R
    buyer_return = base_return

    # Seller's return is -R
    seller_return = -base_return

    # Print the full equation for the buyer's return
    print(f"Buyer's 15s Return = (p' / p) - 1")
    # Using round to avoid floating point representation issues in the output string
    print(f"Calculation: ({p_prime} / {p}) - 1 = {round(buyer_return, 5)}")
    print("-" * 30)

    # Print the full equation for the seller's return
    print(f"Seller's 15s Return = -((p' / p) - 1)")
    print(f"Calculation: -(({p_prime} / {p}) - 1) = {round(seller_return, 5)}")
    print("-" * 30)

    # Verify the condition from the problem description
    is_higher = buyer_return > seller_return
    print(f"\nIs the buyer's return ({round(buyer_return, 5)}) higher than the seller's return ({round(seller_return, 5)})?")
    print(f"Result: {is_higher}")
    print("\nThis confirms the initial observation, which occurs because the subsequent price p' is greater than the trade price p.")

# Example values reflecting the scenario (p' > p)
# Let's assume a trade happened at 100.00
example_trade_price = 100.00
# And 15 seconds later, the mid-price rose to 100.02
example_subsequent_mid_price = 100.02

calculate_returns(example_trade_price, example_subsequent_mid_price)