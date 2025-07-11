import sys

# This script demonstrates why the 15-second return is structurally higher for
# buy-side orders than sell-side orders, based on the provided definitions.
# This phenomenon stems from the mathematical nature of the return calculation
# combined with the bid-ask spread.

def demonstrate_return_difference():
    """
    Calculates and explains the difference in returns for buyers and sellers
    using a sample bid-ask spread.
    """
    # --- Step 1: Define a sample market state ---
    # We set an example bid and ask price for a stock.
    bid_price = 100.00
    ask_price = 100.20
    
    # The bid-ask spread is the difference between ask and bid.
    spread = ask_price - bid_price

    # The price `p'` is the mid-price 15 seconds after the trade.
    # To isolate the structural effect, we assume the mid-price does not change.
    # p_prime is the future mid-price.
    p_prime = (ask_price + bid_price) / 2

    print(f"--- Market State ---")
    print(f"Bid Price: {bid_price:.2f}")
    print(f"Ask Price: {ask_price:.2f}")
    print(f"Spread (Ask - Bid): {spread:.2f}")
    print(f"Future Mid-Price (p'): {p_prime:.3f}\n")
    
    # --- Step 2: Define trade prices and calculate returns ---
    # A buyer-initiated trade executes at the ask price.
    p_buy = ask_price
    # A seller-initiated trade executes at the bid price.
    p_sell = bid_price
    
    # Calculate returns based on the definitions provided.
    # Buyer's return = p'/p - 1
    return_buy = (p_prime / p_buy) - 1
    
    # Seller's return = - (p'/p - 1)
    return_sell = -((p_prime / p_sell) - 1)
    
    print(f"--- Return Calculation ---")
    print(f"Buyer's trade price (p_buy): {p_buy:.2f}")
    print(f"Buyer's return = ({p_prime:.3f} / {p_buy:.2f}) - 1 = {return_buy:.8f}")
    print(f"\nSeller's trade price (p_sell): {p_sell:.2f}")
    print(f"Seller's return = - (({p_prime:.3f} / {p_sell:.2f}) - 1) = {return_sell:.8f}\n")
    
    # --- Step 3: Compare returns and show the mathematical proof ---
    # We demonstrate that the buyer's return is higher than the seller's.
    difference_in_returns = return_buy - return_sell
    
    print(f"--- Comparison ---")
    print(f"Is Buyer's Return > Seller's Return? {'Yes' if return_buy > return_sell else 'No'}")
    print(f"Difference (R_buy - R_sell): {difference_in_returns:.8f}\n")
    
    # The algebraic simplification shows the difference is (ask-bid)^2 / (2*ask*bid)
    # Let's calculate this value to verify our algebra.
    # This equation represents the final equation for the difference.
    numerator = (ask_price - bid_price)**2
    denominator = 2 * ask_price * bid_price
    theoretical_difference = numerator / denominator

    print(f"--- Mathematical Verification ---")
    print("The difference in returns can be shown to be: (ask - bid)^2 / (2 * ask * bid)")
    print(f"Numerator (ask - bid)^2 = ({ask_price:.2f} - {bid_price:.2f})^2 = {numerator:.4f}")
    print(f"Denominator (2 * ask * bid) = 2 * {ask_price:.2f} * {bid_price:.2f} = {denominator:.2f}")
    print(f"Theoretical Difference = {numerator:.4f} / {denominator:.2f} = {theoretical_difference:.8f}\n")

    print(f"Conclusion: The calculated difference ({difference_in_returns:.8f}) matches the theoretical one.")
    print("This shows the buyer's return is structurally higher due to the return definitions and the bid-ask spread.")

if __name__ == "__main__":
    demonstrate_return_difference()
<<<A>>>