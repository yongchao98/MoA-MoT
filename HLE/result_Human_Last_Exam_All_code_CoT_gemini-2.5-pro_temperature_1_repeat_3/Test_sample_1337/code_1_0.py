import random

def solve():
    """
    This function simulates the airline ticket selling process to find the
    booking policy that maximizes expected revenue.
    """

    # --- Problem Parameters ---
    SEAT_CAPACITY = 10
    TICKET_PRICE_CHEAP = 100
    TICKET_PRICE_EXPENSIVE = 200
    NUM_SIMULATIONS = 500000  # Use a large number for accurate results

    # --- Helper function to get customer demand ---
    def get_demand():
        """
        Returns a demand of 0, 1, or 2 based on the probabilities:
        P(0) = 1/4, P(1) = 1/2, P(2) = 1/4
        """
        # random.choices is an efficient way to handle weighted random sampling
        return random.choices([0, 1, 2], weights=[0.25, 0.5, 0.25], k=1)[0]

    # --- Main Simulation Logic ---
    max_expected_revenue = -1.0
    optimal_booking_limit = -1
    optimal_avg_cheap_sold = 0
    optimal_avg_expensive_sold = 0

    print("Running simulation to find the optimal booking policy...")
    # Iterate through all possible booking limits for cheap tickets (b from 0 to 10)
    for b in range(SEAT_CAPACITY + 1):
        total_cheap_sold_for_b = 0
        total_exp_sold_for_b = 0

        # Run the simulation for the current booking limit 'b'
        for _ in range(NUM_SIMULATIONS):
            cheap_sold = 0
            exp_sold = 0

            # Loop through the 14-day selling period
            for day in range(1, 15):
                # Determine daily customer demand
                demand_c1 = get_demand()
                
                # Class 2 customers only appear in the second week (days 8-14)
                demand_c2 = 0
                if day > 7:
                    demand_c2 = get_demand()

                # Process Class 2 customers first (they have priority)
                for _ in range(demand_c2):
                    if (cheap_sold + exp_sold) >= SEAT_CAPACITY:
                        break  # Stop if flight is full
                    
                    if cheap_sold < b:
                        # Sell a cheap ticket if within the booking limit
                        cheap_sold += 1
                    else:
                        # If cheap limit is reached, Class 2 considers an expensive ticket
                        if random.random() < 0.5:
                            exp_sold += 1
                
                # Process Class 1 customers
                for _ in range(demand_c1):
                    if (cheap_sold + exp_sold) >= SEAT_CAPACITY:
                        break  # Stop if flight is full
                    
                    # Class 1 only buys cheap tickets, if available within the limit
                    if cheap_sold < b:
                        cheap_sold += 1

            # Accumulate the number of tickets sold in this simulation run
            total_cheap_sold_for_b += cheap_sold
            total_exp_sold_for_b += exp_sold
        
        # Calculate the average number of tickets sold for this policy
        avg_cheap = total_cheap_sold_for_b / NUM_SIMULATIONS
        avg_exp = total_exp_sold_for_b / NUM_SIMULATIONS
        
        # Calculate the expected revenue for this policy
        current_expected_revenue = (avg_cheap * TICKET_PRICE_CHEAP) + (avg_exp * TICKET_PRICE_EXPENSIVE)

        # Check if this policy is the best one found so far
        if current_expected_revenue > max_expected_revenue:
            max_expected_revenue = current_expected_revenue
            optimal_booking_limit = b
            optimal_avg_cheap_sold = avg_cheap
            optimal_avg_expensive_sold = avg_exp

    # --- Final Output ---
    print("\nSimulation complete. The maximum expected revenue is found.")
    print(f"The optimal policy is to set a booking limit of {optimal_booking_limit} for cheap tickets.")
    print("\nThe final equation for the maximum expected revenue is:")
    print(f"{optimal_avg_cheap_sold:.4f} * {TICKET_PRICE_CHEAP} + {optimal_avg_expensive_sold:.4f} * {TICKET_PRICE_EXPENSIVE} = {max_expected_revenue:.2f}")

    # Return the final numeric answer for the platform
    return max_expected_revenue

# Execute the simulation and print the results
max_revenue_value = solve()
print(f"<<<{max_revenue_value:.2f}>>>")