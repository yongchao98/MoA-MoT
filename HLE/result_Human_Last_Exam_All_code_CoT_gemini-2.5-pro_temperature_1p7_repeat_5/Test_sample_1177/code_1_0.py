def analyze_risk_of_loss():
    """
    This function analyzes the scenario based on the Sale of Goods Act
    to determine who bears the risk of loss.
    """

    # --- Factual timeline from the problem ---
    contract_date = 2 # June 2nd, contract made for the MacBook Pro.
    notification_date = 5 # June 5th, Jake notified Luke the repair was done.
    destruction_date = 6 # Overnight before June 6th, the laptop was destroyed.
    purchase_price = 1000 # The agreed upon price in dollars.

    # --- Legal conditions for transfer of risk ---
    # According to the Sale of Goods Act, for specific goods needing repairs,
    # risk transfers when two conditions are met.

    # Condition 1: Were the goods put into a deliverable state?
    # Jake replaced the screen as promised before the flood.
    seller_completed_repairs = True

    # Condition 2: Did the buyer receive notice that the work was done?
    # Jake texted Luke on June 5th that the laptop was ready.
    buyer_received_notice = True

    # The "equation" to determine if risk passed is checking if both conditions are true.
    # Risk Pass = (Condition 1) AND (Condition 2)
    risk_passed_to_buyer = seller_completed_repairs and buyer_received_notice

    print("Analyzing the transfer of risk for the MacBook Pro transaction:")
    print(f" - On June {contract_date}, a contract was made for ${purchase_price}, but the item required repairs.")
    print(f" - By June {notification_date}, Jake had completed the repairs.")
    print(f" - On June {notification_date}, Jake notified Luke the laptop was ready for pickup.")
    print(f" - Before the pickup on June {destruction_date}, the item was destroyed.")
    print("\n--- Legal Analysis ---")

    if risk_passed_to_buyer:
        print("Conclusion: The risk of loss passed from the seller (Jake) to the buyer (Luke).")
        print("This occurred because both legal conditions were met before the flood:")
        print(f"1. Repairs Completed: {seller_completed_repairs}")
        print(f"2. Buyer Notified: {buyer_received_notice}")
        print("\nTherefore, Jake is not required to return the money. The loss falls on Luke.")
        print("\n--- Correct Answer Explanation ---")
        print("Answer B accurately explains this: Jake is not required to return the money paid to him by Luke because risk of loss passed from Jake to Luke as a result of Jake completing the repairs, and thereby putting the laptop in a deliverable state, and notifying Luke of a precise date on which the laptop will be ready.")
    else:
        # This path is not taken based on the facts.
        print("Conclusion: The risk of loss remained with the seller (Jake).")


# Execute the analysis
analyze_risk_of_loss()
