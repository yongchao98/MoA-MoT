import pandas as pd

def simulate_market_returns():
    """
    Simulates trades in a simple market to demonstrate why a buyer's 15s return
    is typically higher than a seller's.
    """
    # --- 1. Define Market Parameters ---
    # Let's assume a stock with a stable mid-price and a constant bid-ask spread.
    mid_price = 10.00  # The "true" value of the stock in our simple model
    spread = 0.02      # The difference between the best ask and best bid prices

    # --- 2. Calculate Bid and Ask Prices ---
    # An aggressive buyer pays the ask price.
    # An aggressive seller receives the bid price.
    p_ask = mid_price + spread / 2
    p_bid = mid_price - spread / 2

    # --- 3. Model the Future Price ---
    # For the 15s return calculation, we need the price 15s after the trade.
    # Let's assume the price simply reverts to the "true" mid-price.
    # This is a common phenomenon known as "bid-ask bounce".
    p_future = mid_price

    # --- 4. Calculate Returns ---
    # Buyer's return is p'/p - 1
    # Seller's return is -(p'/p - 1)

    # For a buyer-initiated trade:
    p_trade_buy = p_ask
    return_buy = (p_future / p_trade_buy) - 1

    # For a seller-initiated trade:
    p_trade_sell = p_bid
    # Note the formula for the seller is -(p'/p - 1)
    return_sell = -((p_future / p_trade_sell) - 1)

    # --- 5. Print the Explanation and Results ---
    print("### Market Simulation ###")
    print(f"Stable Mid-Price (M): ${mid_price:.2f}")
    print(f"Bid-Ask Spread (S): ${spread:.2f}")
    print(f"--> Buyer Pays (Ask Price): ${p_ask:.2f}")
    print(f"--> Seller Receives (Bid Price): ${p_bid:.2f}")
    print(f"Price 15s later (p'): ${p_future:.2f} (assumed mean-reversion)\n")

    print("### Buyer's Return Calculation ###")
    print(f"Formula: r_buy = p' / p_ask - 1")
    print(f"Equation: r_buy = {p_future:.2f} / {p_ask:.2f} - 1")
    print(f"Result: r_buy = {return_buy:.6f}\n")

    print("### Seller's Return Calculation ###")
    print(f"Formula: r_sell = -(p' / p_bid - 1)")
    print(f"Equation: r_sell = -({p_future:.2f} / {p_bid:.2f} - 1)")
    print(f"Result: r_sell = {return_sell:.6f}\n")

    print("### Comparison ###")
    print(f"Is Buyer's Return > Seller's Return?")
    is_greater = return_buy > return_sell
    print(f"Is {return_buy:.6f} > {return_sell:.6f}?  -->  {is_greater}")
    print("\nConclusion: The phenomenon holds true in this model.")
    print("The core reason is that the buyer's and seller's returns (which represent the cost of crossing the spread) are calculated against different bases (ask price vs. bid price). This is inherent to the return definition.")

simulate_market_returns()
<<<A>>>