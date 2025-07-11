import numpy as np

def market_return_simulation():
    """
    Simulates short-term returns to demonstrate the effect of the bid-ask spread.
    This illustrates why, on average, the post-trade mid-price tends to be
    higher relative to the trade price.
    """
    # 1. Define market parameters for a stock.
    # Let's assume a mid-price of $100 for the stock.
    mid_price = 100.0
    # The half-spread (s) is the difference between the mid-price and the bid/ask.
    # Let's set a half-spread of 5 cents.
    s = 0.05

    # 2. Determine trade and post-trade prices based on the bid-ask spread.
    # A buyer-initiated (aggressive) trade occurs at the ask price.
    p_buy = mid_price + s
    # A seller-initiated (aggressive) trade occurs at the bid price.
    p_sell = mid_price - s

    # For this simple model (called the 'bid-ask bounce'), we assume the price
    # reverts to the original mid-price 15 seconds after the trade.
    p_prime = mid_price

    # 3. Calculate the 15s return for each type of trade.
    # Formula: p'/p - 1
    return_for_buy_trade = (p_prime / p_buy) - 1
    return_for_sell_trade = (p_prime / p_sell) - 1

    # 4. Calculate the average return, assuming an equal number of buy- and sell-initiated trades.
    # The average is simply the average of the two return types.
    average_return = (return_for_buy_trade + return_for_sell_trade) / 2

    # 5. Print the explanation and the final equation with numbers.
    print("This simulation shows why the average 15s return is typically positive.")
    print("The reason is the market's bid-ask spread structure.\n")
    print("--- Model Parameters ---")
    print(f"Assumed Mid-Price (p' after bounce): {p_prime}")
    print(f"Half-Spread (s): {s}")
    print(f"Buy Trade Price (p_buy = mid + s): {p_buy}")
    print(f"Sell Trade Price (p_sell = mid - s): {p_sell}\n")

    print("--- Individual Return Calculations ---")
    print(f"Return for a Buy Trade ({p_prime}/{p_buy} - 1): {return_for_buy_trade:.8f}")
    print(f"Return for a Sell Trade ({p_prime}/{p_sell} - 1): {return_for_sell_trade:.8f}\n")
    print("Note that the positive return for a sell trade is larger in magnitude than the negative return for a buy trade.\n")
    
    print("--- Final Average Return Equation ---")
    # We output the final equation showing all the numbers, as requested.
    equation_string = (
        f"Average Return = (Return for Buy + Return for Sell) / 2\n"
        f"Average Return = ({return_for_buy_trade:.8f} + {return_for_sell_trade:.8f}) / 2 = {average_return:.8f}"
    )
    print(equation_string)
    
    print("\nThe simulation yields a positive average return. This phenomenon is a feature of all limit-order markets.")
    print("The most fundamental economic cause for this market structure is Information Asymmetry.")

# Run the simulation.
market_return_simulation()
<<<C>>>