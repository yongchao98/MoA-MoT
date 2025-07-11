import random

def simulate_market_returns():
    """
    This simulation demonstrates how asymmetric price impact, stemming from different
    motivations of buyers and sellers, leads to the observed return pattern.
    """
    # --- Simulation Parameters ---
    # We model the hypothesis that buyers are often more "information-driven" (larger price impact)
    # than sellers, who may be "liquidity-driven" (smaller price impact).
    p_mid = 100.0  # Initial mid-price of the stock
    spread = 0.02  # The difference between the best ask and bid prices
    buy_impact = 0.05  # The average price increase 15s after a buy-initiated trade
    sell_impact = -0.03 # The average price decrease 15s after a sell-initiated trade (note: abs(sell_impact) < buy_impact)
    num_trades = 10000 # Number of trades to simulate for a stable average

    # --- Data storage ---
    buy_side_order_returns = []
    sell_side_order_returns = []

    # --- Simulation Loop ---
    for _ in range(num_trades):
        # We randomly simulate a trade being initiated by either a buyer or a seller
        if random.random() > 0.5:
            # Case 1: A "buy-side order" (an aggressive buyer executing a trade)
            # The trade price 'p' is the ask price
            p = p_mid + spread / 2
            
            # The future price 'p_prime' reflects the strong positive impact of the buy
            p_prime = p_mid + buy_impact
            
            # Per definition, the 15s return for the buyer is p'/p - 1
            return_for_buyer = (p_prime / p) - 1
            buy_side_order_returns.append(return_for_buyer)
        else:
            # Case 2: A "sell-side order" (an aggressive seller executing a trade)
            # The trade price 'p' is the bid price
            p = p_mid - spread / 2
            
            # The future price 'p_prime' reflects the weaker negative impact of the sell
            p_prime = p_mid + sell_impact
            
            # Per definition, the 15s return for the seller is -(p'/p - 1)
            return_for_seller = -((p_prime / p) - 1)
            sell_side_order_returns.append(return_for_seller)

    # --- Calculate and Report Final Averages ---
    avg_buy_side_return = sum(buy_side_order_returns) / len(buy_side_order_returns)
    avg_sell_side_return = sum(sell_side_order_returns) / len(sell_side_order_returns)

    print("This simulation models returns based on asymmetric price impact.")
    print(f"Key Assumptions: buy_impact = {buy_impact}, sell_impact = {sell_impact}")
    print("-" * 45)
    print("Final Calculation of Average Returns:")
    print(f"Average Return for Buy-Side Orders = {avg_buy_side_return:.6f}")
    print(f"Average Return for Sell-Side Orders = {avg_sell_side_return:.6f}")
    print("-" * 45)
    
    if avg_buy_side_return > avg_sell_side_return:
        print("Conclusion: The simulation confirms that the average return for buy-side orders is higher.")
    else:
        print("Conclusion: The simulation does not show higher returns for buy-side orders.")

simulate_market_returns()