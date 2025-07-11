def simulate_asymmetric_returns():
    """
    This simulation demonstrates how different motivations of buyers and sellers,
    leading to asymmetric price impacts, can result in the 15s return for
    buy-side orders being higher than for sell-side orders.
    """
    # --- Step 1: Set up market parameters ---
    # Assume a stock is trading around $100
    p_bid = 99.95  # The highest price someone is willing to pay
    p_ask = 100.05 # The lowest price someone is willing to sell for
    print(f"Initial State: Bid Price = {p_bid}, Ask Price = {p_ask}\n")

    # --- Step 2: Define asymmetric price impact based on different motivations ---
    # We assume buy orders have a larger price impact than sell orders.
    # This reflects that buyers might be more consistently driven by positive information,
    # while sellers have mixed motivations (negative info, liquidity needs, etc.).
    delta_up = 0.10   # Price moves up by 10 cents after an aggressive buy
    delta_down = 0.06 # Price moves down by 6 cents after an aggressive sell
    print(f"Assumption: Upward price impact (delta_up) = {delta_up}")
    print(f"Assumption: Downward price impact (delta_down) = {delta_down}")
    print("Note: delta_up > delta_down, modeling asymmetric impact.\n")

    # --- Step 3: Calculate return for a Buyer-Initiated Trade ---
    print("--- Calculating Return for Buyer-Initiated Order ---")
    p_buy_trade = p_ask
    print(f"A buyer initiates a trade. The trade price 'p' is the ask price: {p_buy_trade}")
    
    # The price 15 seconds later reflects the upward impact
    p_prime_after_buy = p_buy_trade + delta_up
    print(f"Price 15s later 'p'' is the trade price + delta_up: {p_buy_trade} + {delta_up} = {p_prime_after_buy}")

    # Calculate the buyer's 15s return
    # Formula: p'/p - 1
    return_buy = p_prime_after_buy / p_buy_trade - 1
    print(f"Buyer's 15s Return = p'/p - 1 = {p_prime_after_buy} / {p_buy_trade} - 1 = {return_buy:.6f}\n")

    # --- Step 4: Calculate return for a Seller-Initiated Trade ---
    print("--- Calculating Return for Seller-Initiated Order ---")
    p_sell_trade = p_bid
    print(f"A seller initiates a trade. The trade price 'p' is the bid price: {p_sell_trade}")

    # The price 15 seconds later reflects the downward impact
    p_prime_after_sell = p_sell_trade - delta_down
    print(f"Price 15s later 'p'' is the trade price - delta_down: {p_sell_trade} - {delta_down} = {p_prime_after_sell}")

    # Calculate the seller's 15s return
    # Formula: -(p'/p - 1) which is 1 - p'/p
    return_sell = 1 - (p_prime_after_sell / p_sell_trade)
    print(f"Seller's 15s Return = 1 - p'/p = 1 - {p_prime_after_sell} / {p_sell_trade} = {return_sell:.6f}\n")

    # --- Step 5: Compare the returns ---
    print("--- Final Comparison ---")
    print(f"Buyer's Return: {return_buy:.6f}")
    print(f"Seller's Return: {return_sell:.6f}")
    if return_buy > return_sell:
        print("Conclusion: The buyer's return is indeed higher than the seller's return in this scenario.")
    else:
        print("Conclusion: The buyer's return is NOT higher than the seller's return.")

# Execute the simulation
simulate_asymmetric_returns()