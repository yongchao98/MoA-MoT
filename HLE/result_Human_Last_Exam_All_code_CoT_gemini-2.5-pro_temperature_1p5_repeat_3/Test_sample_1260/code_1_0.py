import numpy as np

def simulate_trading_returns():
    """
    Simulates trading and calculates returns to demonstrate why buy-side
    returns are often higher than sell-side returns.
    """
    # --- Market Parameters ---
    initial_mid_price = 100.0
    spread = 0.10
    ask_price = initial_mid_price + spread / 2
    bid_price = initial_mid_price - spread / 2
    num_trades = 2000

    # --- Asymmetric Price Impact Assumption (The "Why") ---
    # We model that buy orders have a larger average price impact,
    # reflecting that they are, on average, more "informed".
    # Price impact is how much the mid-price moves 15s after a trade.
    avg_buy_impact = 0.03   # Price moves up by 0.03 on avg after a buy
    avg_sell_impact = -0.02  # Price moves down by 0.02 on avg after a sell
    impact_std_dev = 0.05    # Volatility of impact

    buy_side_returns = []
    sell_side_returns = []
    
    print("Simulating trades based on the principle of asymmetric price impact...\n")

    for i in range(num_trades):
        # Randomly choose whether it's a buy or sell initiated trade
        if np.random.rand() > 0.5:
            # --- Buy-initiated Trade ---
            trade_price = ask_price
            # The price moves up due to the buy order's impact
            price_impact = np.random.normal(avg_buy_impact, impact_std_dev)
            price_15s_later = initial_mid_price + price_impact
            
            # Calculate buyer's return
            buyer_return = price_15s_later / trade_price - 1
            buy_side_returns.append(buyer_return)

            # For the first buy trade, print the equation
            if len(buy_side_returns) == 1:
                print("--- Sample Calculation for a Buy-Side Order ---")
                print(f"Trade at Ask Price (p): {trade_price:.4f}")
                print(f"Mid-Price 15s later (p'): {price_15s_later:.4f}")
                print(f"Buyer's Return Formula: p'/p - 1")
                print(f"Calculation: {price_15s_later:.4f} / {trade_price:.4f} - 1 = {buyer_return:.6f}\n")

        else:
            # --- Sell-initiated Trade ---
            trade_price = bid_price
            # The price moves down due to the sell order's impact
            price_impact = np.random.normal(avg_sell_impact, impact_std_dev)
            price_15s_later = initial_mid_price + price_impact
            
            # Calculate seller's return
            seller_return = -(price_15s_later / trade_price - 1)
            sell_side_returns.append(seller_return)
            
            # For the first sell trade, print the equation
            if len(sell_side_returns) == 1:
                print("--- Sample Calculation for a Sell-Side Order ---")
                print(f"Trade at Bid Price (p): {trade_price:.4f}")
                print(f"Mid-Price 15s later (p'): {price_15s_later:.4f}")
                print(f"Seller's Return Formula: -(p'/p - 1)")
                print(f"Calculation: -({price_15s_later:.4f} / {trade_price:.4f} - 1) = {seller_return:.6f}\n")

    # --- Calculate and Print Final Averages ---
    avg_buy_return = np.mean(buy_side_returns)
    avg_sell_return = np.mean(sell_side_returns)

    print("--- Final Results ---")
    print(f"Number of buy-side trades: {len(buy_side_returns)}")
    print(f"Number of sell-side trades: {len(sell_side_returns)}")
    print(f"Average 15s return for buy-side orders: {avg_buy_return:.6f}")
    print(f"Average 15s return for sell-side orders: {avg_sell_return:.6f}")
    
    if avg_buy_return > avg_sell_return:
        print("\nConclusion: The simulation confirms the observation that average buy-side returns are higher.")
    else:
        print("\nConclusion: The simulation does not confirm the observation (this may happen by chance).")

# Run the simulation
simulate_trading_returns()
