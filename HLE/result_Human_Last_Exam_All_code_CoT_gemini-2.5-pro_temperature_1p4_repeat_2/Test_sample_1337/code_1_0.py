import random

def solve_flight_revenue():
    """
    Finds the maximum expected revenue for the flight by simulating different
    booking policies based on a protection level.
    """
    
    def get_demand():
        """Returns customer demand (0, 1, or 2) based on the given probability distribution."""
        rand_num = random.random()
        if rand_num < 0.25:
            return 0
        elif rand_num < 0.75:  # 0.25 (for 0) + 0.50 (for 1)
            return 1
        else:
            return 2

    def simulate_one_booking_period(protection_level):
        """
        Simulates one 14-day booking period for a given protection level.
        A protection level is the number of seats reserved for expensive fares.
        
        Returns a tuple of: (total_revenue, num_cheap_sold, num_expensive_sold)
        """
        capacity = 10
        price_cheap = 100
        price_expensive = 200
        limit_cheap_tickets = capacity - protection_level

        sold_cheap = 0
        sold_expensive = 0
        
        # Simulate each of the 14 days
        for day in range(1, 15):
            # Determine demand for the day
            c1_demand = get_demand()
            c2_demand = get_demand() if day > 7 else 0

            # Process Class 2 customers first (they have priority for cheap tickets)
            for _ in range(c2_demand):
                if sold_cheap + sold_expensive >= capacity:
                    break  # Flight is full
                
                if sold_cheap < limit_cheap_tickets:
                    # Sell a cheap ticket
                    sold_cheap += 1
                else:
                    # No cheap tickets left, 50% chance to buy an expensive one
                    if random.random() < 0.5:
                        sold_expensive += 1

            # Process Class 1 customers
            for _ in range(c1_demand):
                if sold_cheap + sold_expensive >= capacity:
                    break  # Flight is full
                
                # Class 1 only buys cheap tickets
                if sold_cheap < limit_cheap_tickets:
                    sold_cheap += 1
        
        revenue = (sold_cheap * price_cheap) + (sold_expensive * price_expensive)
        return revenue, sold_cheap, sold_expensive

    # --- Main execution ---
    num_simulations = 200000
    policy_stats = {}  # To store results for each policy k

    # Iterate through all possible policies (protection levels k from 0 to 10)
    for k in range(11):
        total_revenue = 0
        total_cheap_sold = 0
        total_expensive_sold = 0
        
        for _ in range(num_simulations):
            revenue, cheap, expensive = simulate_one_booking_period(k)
            total_revenue += revenue
            total_cheap_sold += cheap
            total_expensive_sold += expensive
        
        # Store the average results for policy k
        policy_stats[k] = {
            'expected_revenue': total_revenue / num_simulations,
            'avg_cheap_sold': total_cheap_sold / num_simulations,
            'avg_expensive_sold': total_expensive_sold / num_simulations
        }

    # Find the policy with the maximum expected revenue
    best_protection_level = -1
    max_expected_revenue = -1
    for k, stats in policy_stats.items():
        if stats['expected_revenue'] > max_expected_revenue:
            max_expected_revenue = stats['expected_revenue']
            best_protection_level = k

    # Get the detailed stats for the best policy
    best_policy_stats = policy_stats[best_protection_level]
    avg_cheap_sold = best_policy_stats['avg_cheap_sold']
    avg_expensive_sold = best_policy_stats['avg_expensive_sold']
    final_revenue = avg_cheap_sold * 100 + avg_expensive_sold * 200

    # Print the explanation and final result
    print(f"To find the maximum revenue, we simulate {num_simulations} selling periods for each possible booking policy.")
    print("A policy is defined by the 'protection level': the number of seats reserved for expensive tickets.")
    print("\n--- Optimal Policy Found ---")
    print(f"The best policy is to protect {best_protection_level} seats for expensive fares.")
    print(f"This means a maximum of {10 - best_protection_level} cheap tickets will be offered for sale.")
    
    print("\n--- Expected Results for the Optimal Policy ---")
    print(f"Average cheap tickets sold: {avg_cheap_sold:.2f}")
    print(f"Average expensive tickets sold: {avg_expensive_sold:.2f}")
    
    print("\n--- Final Revenue Calculation ---")
    print(f"The final maximum expected revenue is calculated from the average sales:")
    print(f"({avg_cheap_sold:.2f} * $100) + ({avg_expensive_sold:.2f} * $200) = ${final_revenue:.2f}")

    print(f"\n<<<{final_revenue:.2f}>>>")

solve_flight_revenue()