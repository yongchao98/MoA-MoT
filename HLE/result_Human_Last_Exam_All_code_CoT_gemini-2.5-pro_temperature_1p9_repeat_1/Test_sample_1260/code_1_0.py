import numpy as np

def market_simulation():
    """
    Simulates trades in a market to demonstrate how different motivations
    of buyers and sellers can lead to the observed return asymmetry.

    The core of the simulation is the assumption that the average price impact of
    a buy-initiated trade is greater than the magnitude of the average price
    impact of a sell-initiated trade.
    """

    # --- Market & Simulation Parameters ---
    p_mid_initial = 100.0
    spread = 0.02
    num_trades = 20000

    # --- Asymmetric Price Impact Assumption ---
    # This reflects that buyers' motivations might lead to stronger price moves.
    # Impact is the change in mid-price 15s after the trade.
    buy_impact_mean = 0.03
    sell_impact_magnitude_mean = 0.02 # Sells have a smaller price impact
    impact_std = 0.01

    buy_initiator_returns = []
    sell_initiator_returns = []

    print("--- Simulation Setup ---")
    print("This simulation models a market where buy orders have a larger price impact than sell orders.")
    print(f"Assumed Average Buy Impact: +{buy_impact_mean}")
    print(f"Assumed Average Sell Impact: -{sell_impact_magnitude_mean}")
    print("--------------------------\n")

    for i in range(num_trades):
        p_ask = p_mid_initial + spread / 2
        p_bid = p_mid_initial - spread / 2

        # In each trial, randomly decide if it is a buy- or sell-initiated trade
        if np.random.rand() > 0.5:
            # Case 1: Buy-initiated trade
            p = p_ask
            impact = np.random.normal(buy_impact_mean, impact_std)
            p_prime = p_mid_initial + impact
            # The initiator is the buyer. Their return is p'/p - 1.
            initiator_return = (p_prime / p) - 1
            buy_initiator_returns.append(initiator_return)
        else:
            # Case 2: Sell-initiated trade
            p = p_bid
            impact_magnitude = np.random.normal(sell_impact_magnitude_mean, impact_std)
            p_prime = p_mid_initial - impact_magnitude
            # The initiator is the seller. Their return is -(p'/p - 1).
            initiator_return = -((p_prime / p) - 1)
            sell_initiator_returns.append(initiator_return)

    # Calculate the average returns for each type of initiator
    avg_return_buy = np.mean(buy_initiator_returns)
    avg_return_sell = np.mean(sell_initiator_returns)

    # The "final equation" from the problem is the comparison of these two average returns.
    # We print each number in the equation.
    print("--- Final Equation Check ---")
    print("Is the 15s return higher for buy-side orders than sell-side orders?")
    print("This corresponds to: Avg(Return for Buy Initiator) > Avg(Return for Sell Initiator)\n")
    print("Calculated Values:")
    print(f"Avg Return for Buy Initiator : {avg_return_buy}")
    print(f"Avg Return for Sell Initiator: {avg_return_sell}\n")

    print("Final Result:")
    result_is_true = avg_return_buy > avg_return_sell
    print(f"The inequality {avg_return_buy} > {avg_return_sell} is {result_is_true}.")

market_simulation()
<<<L>>>