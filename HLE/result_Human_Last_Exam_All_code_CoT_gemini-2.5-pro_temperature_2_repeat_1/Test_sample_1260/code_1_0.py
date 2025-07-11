def explain_return_asymmetry():
    """
    Analyzes and explains the asymmetry in 15s returns for buyers and sellers
    based on the structure of the market and the definition of returns.
    """
    print("Analyzing the stock return phenomenon step-by-step.")
    print("--------------------------------------------------\n")

    # Step 1: Define a hypothetical market state (pre-trade)
    bid_price = 10.00
    ask_price = 10.02
    pre_trade_mid_price = (bid_price + ask_price) / 2
    
    print(f"Let's assume a stock has a pre-trade Bid Price of ${bid_price:.2f} and an Ask Price of ${ask_price:.2f}.")
    print(f"The pre-trade mid-price is (${bid_price:.2f} + ${ask_price:.2f}) / 2 = ${pre_trade_mid_price:.3f}.\n")

    # For this demonstration, we assume the price 15s later (p') reverts to the original mid-price.
    # This helps isolate the structural effect from market trends like momentum.
    post_trade_mid_price = pre_trade_mid_price
    print(f"To find the fundamental reason, let's assume a neutral market where the price")
    print(f"15 seconds after a trade (p') reverts to the original mid-price: ${post_trade_mid_price:.3f}.\n")

    # Step 2: Analyze a Buy-Side Order (where the buyer is the initiator)
    print("--- Case 1: A Buy-Side Initiated Order ---")
    p_buy = ask_price
    print(f"An aggressive buyer crosses the spread. The trade price (p) is the ask price: ${p_buy:.2f}.")
    
    # Calculate buyer's 15s return
    buyer_return = (post_trade_mid_price / p_buy) - 1
    
    print("The buyer's 15s return is calculated as: p' / p - 1")
    print(f"Equation: (${post_trade_mid_price:.3f} / ${p_buy:.2f}) - 1 = {buyer_return:.6f}")
    print(f"This small negative return represents the buyer's cost for getting the trade done immediately.\n")

    # Step 3: Analyze a Sell-Side Order (where the seller is the initiator)
    print("--- Case 2: A Sell-Side Initiated Order ---")
    p_sell = bid_price
    print(f"An aggressive seller crosses the spread. The trade price (p) is the bid price: ${p_sell:.2f}.")
    
    # Calculate seller's 15s return
    seller_return = -((post_trade_mid_price / p_sell) - 1)
    
    print("The seller's 15s return is calculated as: -(p' / p - 1)")
    print(f"Equation: -((${post_trade_mid_price:.3f} / ${p_sell:.2f}) - 1) = {seller_return:.6f}")
    print(f"This small negative return represents the seller's cost for liquidity.\n")

    # Step 4: Compare the two returns
    print("--- Comparison ---")
    print(f"Return for the buy-side order initiator: {buyer_return:.6f}")
    print(f"Return for the sell-side order initiator: {seller_return:.6f}")
    
    is_higher = buyer_return > seller_return
    print(f"\nIs the buyer's return higher than the seller's return? {is_higher}.")
    print(f"The calculation shows {buyer_return:.6f} > {seller_return:.6f}.")

    print("\n--- Conclusion ---")
    print("The buyer's return is higher than the seller's return due to the calculation's structure.")
    print("A buyer's return is benchmarked against the higher ask price, while a seller's is benchmarked against the lower bid price.")
    print("This is a structural property of any market with a bid-ask spread and is not dependent on other factors like momentum, although they can coexist.")
    print("Therefore, this phenomenon is expected in all listed markets (Shanghai, Nasdaq, HK), and the most probable reason is the nature of the 15s return definition itself.")

# Execute the analysis
explain_return_asymmetry()