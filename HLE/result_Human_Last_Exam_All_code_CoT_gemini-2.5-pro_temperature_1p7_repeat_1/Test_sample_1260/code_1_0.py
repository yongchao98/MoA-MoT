def explain_return_asymmetry():
    """
    Illustrates why the average 15s return can be positive.
    The key is an asymmetry between the impact of buy- and sell-initiated trades.
    """
    print("Analyzing the 15s return calculation...")
    print("Let R = p'/p - 1.")
    print("Buyer's Return = R")
    print("Seller's Return = -R")
    print("The observation 'Average Buyer's Return > Average Seller's Return' implies Average(R) > 0.")
    print("\nThis happens if there is an asymmetry. Let's model an example where buy trades have a stronger price impact.")

    # --- Example Calculation with Asymmetric Impact ---
    
    # Trade 1: A buy-initiated trade
    p_buy = 100.20
    # The trade causes a price impact, let's say the mid-price drifts up by 0.15
    p_prime_buy = p_buy + 0.15 
    # Let's write p_prime_buy with the correct number
    p_prime_buy = 100.35
    R_buy = p_prime_buy / p_buy - 1
    
    print("\n--- Trade 1 (Buy-initiated) ---")
    print(f"Initial Trade Price (p): {p_buy}")
    print(f"Mid-price after 15s (p'): {p_prime_buy}")
    print(f"Return Term (p'/p - 1) = {p_prime_buy}/{p_buy} - 1 = {R_buy:.6f}")

    # Trade 2: A sell-initiated trade
    p_sell = 100.00
    # The sell trade causes a smaller downward drift, say -0.10
    p_prime_sell = p_sell - 0.10
    # Let's write p_prime_sell with the correct number
    p_prime_sell = 99.90
    R_sell = p_prime_sell / p_sell - 1

    print("\n--- Trade 2 (Sell-initiated) ---")
    print(f"Initial Trade Price (p): {p_sell}")
    print(f"Mid-price after 15s (p'): {p_prime_sell}")
    print(f"Return Term (p'/p - 1) = {p_prime_sell}/{p_sell} - 1 = {R_sell:.6f}")
    
    # Calculate the average return term
    average_R = (R_buy + R_sell) / 2
    
    print("\n--- Overall Average ---")
    print(f"Average Return Term (Average(R)) = ({R_buy:.6f} + {R_sell:.6f}) / 2 = {average_R:.6f}")

    if average_R > 0:
        print("\nBecause the Average(R) is positive, the Average Buyer's Return is positive, and higher than the Average Seller's Return.")
    else:
        print("\nIn this example, the Average(R) is not positive. The asymmetry would need to be different for the observation to hold.")

# Run the explanation
explain_return_asymmetry()