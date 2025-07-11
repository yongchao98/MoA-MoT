def analyze_15s_return():
    """
    Analyzes the 15s return phenomenon based on the bid-ask spread.

    This function demonstrates that even with no price momentum (i.e., the price on average
    returns to the midpoint), the 15s return for a buyer is, on average, higher than
    for a seller. This is a mathematical consequence of the bid-ask spread.
    """

    # --- Setup with example market data ---
    # Let's assume a stock is trading around $100
    mid_price = 100.00
    # Let's assume a tight spread of 2 cents
    spread = 0.02
    half_spread = spread / 2

    # Calculate Bid and Ask prices
    ask_price = mid_price + half_spread  # Price for an aggressive BUYER
    bid_price = mid_price - half_spread  # Price for an aggressive SELLER

    # --- Assumption for p' ---
    # Let's assume after 15 seconds, the price has no specific trend (no momentum)
    # and on average reverts to the original mid-price.
    p_prime = mid_price

    print("--- Market Assumptions ---")
    print(f"Mid Price: {mid_price:.4f}")
    print(f"Spread: {spread:.4f}")
    print(f"Ask Price (p for aggressive buy): {ask_price:.4f}")
    print(f"Bid Price (p for aggressive sell): {bid_price:.4f}")
    print(f"Assumed Future Mid-Price (p'): {p_prime:.4f}\n")

    # --- Calculate Returns for each trade type ---
    print("--- Return Calculation for Different Scenarios ---")

    # Scenario 1: An aggressive buyer initiates the trade.
    # The trade happens at the Ask price. This person is the "buyer".
    p_buy_initiated = ask_price
    return_for_buyer_in_buy_trade = (p_prime / p_buy_initiated) - 1
    return_for_seller_in_buy_trade = -return_for_buyer_in_buy_trade

    print(f"For a BUYER-initiated trade at p = {p_buy_initiated:.4f}:")
    print(f"  - Buyer's 15s Return (p'/p - 1): {return_for_buyer_in_buy_trade:+.8f}")
    print(f"  - Seller's 15s Return -(p'/p - 1): {return_for_seller_in_buy_trade:+.8f}")
    print("  Here, the seller's return is higher.\n")

    # Scenario 2: An aggressive seller initiates the trade.
    # The trade happens at the Bid price. The counterparty is the "buyer".
    p_sell_initiated = bid_price
    return_for_buyer_in_sell_trade = (p_prime / p_sell_initiated) - 1
    return_for_seller_in_sell_trade = -return_for_buyer_in_sell_trade

    print(f"For a SELLER-initiated trade at p = {p_sell_initiated:.4f}:")
    print(f"  - Buyer's 15s Return (p'/p - 1): {return_for_buyer_in_sell_trade:+.8f}")
    print(f"  - Seller's 15s Return -(p'/p - 1): {return_for_seller_in_sell_trade:+.8f}")
    print("  Here, the buyer's return is higher.\n")

    # --- Calculate the Average Return ---
    # To find the "usual" or average outcome, let's assume aggressive buys and sells
    # are equally likely (50/50 mix). We are interested in the average return for
    # the buyer's side of any trade.
    print("--- Average Return Calculation (assuming 50/50 mix of trades) ---")
    # The buyer's return is either from a buy-initiated trade or a sell-initiated trade
    avg_buyer_return = 0.5 * return_for_buyer_in_buy_trade + 0.5 * return_for_buyer_in_sell_trade
    avg_seller_return = 0.5 * return_for_seller_in_buy_trade + 0.5 * return_for_seller_in_sell_trade

    print(f"Average Buyer's Return: {avg_buyer_return:+.8f}")
    print(f"Average Seller's Return: {avg_seller_return:+.8f}")

    # --- Final Conclusion ---
    is_higher = avg_buyer_return > avg_seller_return
    print("\n--- Conclusion ---")
    print(f"Is the average buyer's return higher than the average seller's return? {is_higher}")
    print("The average buyer's return is positive, while the seller's is negative.")
    print("This happens due to the math of averaging returns across the bid-ask spread.")
    print("This is a structural feature of the market and the return definition, not necessarily due to momentum or sentiment.")
    
    # We can also prove this with the derived formula: s^2 / (Mid^2 - s^2)
    # This formula gives the average of (p'/p - 1)
    s = half_spread
    M = mid_price
    theoretical_avg_return = (s**2) / (M**2 - s**2)
    print(f"\nTheoretical average of (p'/p - 1) is s^2 / (M^2 - s^2) = {theoretical_avg_return:+.8f}")
    print("This value is positive, confirming the observation.")


analyze_15s_return()