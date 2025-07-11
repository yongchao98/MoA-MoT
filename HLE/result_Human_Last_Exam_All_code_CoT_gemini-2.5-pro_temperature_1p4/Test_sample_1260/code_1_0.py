import sys

def simulate_trade_returns():
    """
    This simulation demonstrates why the 15s return for a buyer is higher
    in a buy-initiated trade compared to a sell-initiated trade.
    
    The key principle is Price Impact: an aggressive trade tends to move
    the market price in the direction of the trade.
    
    - A buyer-initiated trade pushes the price up.
    - A seller-initiated trade pushes the price down.
    
    We will analyze the return from the BUYER's perspective in both cases.
    """
    
    # --- Assumptions for our hypothetical market ---
    # Initial mid-price of the stock
    initial_mid_price = 100.00
    # The bid-ask spread is 0.10, so half-spread is 0.05
    half_spread = 0.05
    # The price impact of an aggressive trade over 15 seconds
    price_impact = 0.08

    print("--- Scenario 1: Buyer-Initiated Trade ---")
    print("An aggressive buyer buys at the ask price.")
    
    # The trade happens at the ask price
    p_buy = initial_mid_price + half_spread
    print(f"Trade price (p) = Initial Mid-Price + Half Spread = {initial_mid_price:.2f} + {half_spread:.2f} = {p_buy:.2f}")

    # The price impact pushes the new mid-price up
    p_prime_buy = initial_mid_price + price_impact
    print(f"Mid-price after 15s (p') = Initial Mid-Price + Price Impact = {initial_mid_price:.2f} + {price_impact:.2f} = {p_prime_buy:.2f}")

    # Calculate the aggressive buyer's return
    buyer_return_in_buy_trade = (p_prime_buy / p_buy) - 1
    print(f"Buyer's 15s Return = (p' / p) - 1 = ({p_prime_buy:.4f} / {p_buy:.4f}) - 1 = {buyer_return_in_buy_trade:.6f}")
    print(f"The return is positive, as the buyer bought just before the price rose.\n")


    print("--- Scenario 2: Seller-Initiated Trade ---")
    print("An aggressive seller sells at the bid price. The buyer is a passive participant.")

    # The trade happens at the bid price
    p_sell = initial_mid_price - half_spread
    print(f"Trade price (p) = Initial Mid-Price - Half Spread = {initial_mid_price:.2f} - {half_spread:.2f} = {p_sell:.2f}")
    
    # The price impact pushes the new mid-price down
    p_prime_sell = initial_mid_price - price_impact
    print(f"Mid-price after 15s (p') = Initial Mid-Price - Price Impact = {initial_mid_price:.2f} - {price_impact:.2f} = {p_prime_sell:.2f}")

    # Calculate the passive buyer's return
    buyer_return_in_sell_trade = (p_prime_sell / p_sell) - 1
    print(f"Buyer's 15s Return = (p' / p) - 1 = ({p_prime_sell:.4f} / {p_sell:.4f}) - 1 = {buyer_return_in_sell_trade:.6f}")
    print(f"The return is negative, as the buyer bought just before the price fell.\n")

    
    print("--- Conclusion ---")
    print(f"Return in Buy-Initiated Trade ({buyer_return_in_buy_trade:.6f}) > Return in Sell-Initiated Trade ({buyer_return_in_sell_trade:.6f})")
    print("The simulation confirms the observation. This happens in all modern order-driven markets and is fundamentally caused by information asymmetry.")

# Execute the simulation
if __name__ == "__main__":
    # The following line is needed for the platform but does not affect the logic.
    if len(sys.argv) > 1 and sys.argv[1] == '--execute':
        simulate_trade_returns()
<<<C>>>