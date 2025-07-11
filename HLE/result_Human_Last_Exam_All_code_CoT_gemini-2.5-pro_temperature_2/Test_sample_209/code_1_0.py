import math

# State Variables
time = 0
goldfish_in_pond = 10
shark_eaten_count = 0
shark_last_meal_time = 0

fisherman_caught_count = 0
goldfish_in_basket = 0
basket_capacity = 4

fisherman_is_traveling = False
travel_return_time = -1

# Simulation constants
CATCH_TIME_ROD_B = 5
SHARK_EAT_INTERVAL = 10
BASKET_TRIP_DURATION = 20
SHARK_SURVIVAL_TIME = 11

print("Running a simulation to find the earliest time all fish are gone.")
print("Optimal Strategy: The fisherman uses the fastest rod (Rod B) continuously.")
print("-" * 60)
print("Timeline of Events:")
print(f"Time {time}: START. Pond has {goldfish_in_pond} goldfish. Fisherman begins using Rod B.")

# The simulation runs until all free-swimming goldfish are gone.
while goldfish_in_pond > 0:
    time += 1
    
    # Event: A fisherman event (catching or returning) can happen.
    # Check if the fisherman is currently fishing.
    if not fisherman_is_traveling:
        # Check if it's time to catch a fish with Rod B.
        if time % CATCH_TIME_ROD_B == 0:
            if goldfish_in_pond > 0:
                goldfish_in_pond -= 1
                fisherman_caught_count += 1
                goldfish_in_basket += 1
                print(f"Time {time}: Fisherman catches a goldfish. Pond: {goldfish_in_pond}, Basket: {goldfish_in_basket}.")
            
            # After catching, check if the basket is now full.
            if goldfish_in_basket == basket_capacity:
                fisherman_is_traveling = True
                travel_return_time = time + BASKET_TRIP_DURATION
                print(f"Time {time}: Basket is full. Fisherman leaves for a new basket (returns at T={travel_return_time}).")
    
    # Check if the fisherman returns from his trip.
    elif time == travel_return_time:
        fisherman_is_traveling = False
        goldfish_in_basket = 0  # He has a new, empty basket.
        print(f"Time {time}: Fisherman returns and resumes fishing with Rod B.")

    # Event: A shark event can happen.
    # Check if it's time for the shark to eat.
    # Note: In this scenario, the shark never eats enough fish to speed up.
    if time > 0 and time % SHARK_EAT_INTERVAL == 0:
        if goldfish_in_pond > 0:
            goldfish_in_pond -= 1
            shark_eaten_count += 1
            shark_last_meal_time = time
            print(f"Time {time}: Shark eats a free-swimming goldfish. Pond: {goldfish_in_pond}.")

# End of loop: All free-swimming goldfish have been removed from the pond.
final_time_all_goldfish_gone = time
print("-" * 60)
print(f"Time {final_time_all_goldfish_gone}: All free-swimming goldfish are now gone from the pond.")
print(f"  - Total goldfish caught by fisherman: {fisherman_caught_count}")
print(f"  - Total goldfish eaten by shark: {shark_eaten_count}")
print(f"\nThe shark's last meal was at T={shark_last_meal_time} minutes.")

# Determine when the shark dies.
# The rule for extended survival time is: "eaten more than 4 fish within the last 60 minutes".
# At T=40 (its last meal), the shark has eaten exactly 4 fish, so the rule does not apply.
print(f"The shark's survival time without food is {SHARK_SURVIVAL_TIME} minutes.")
print("Since there are no goldfish left, the shark will starve.")

final_time = shark_last_meal_time + SHARK_SURVIVAL_TIME

print("\nFinal Calculation:")
print("The time when there are no more fish left is determined by when the shark starves.")
print(f"Earliest Time = (Shark's Last Meal Time) + (Shark's Survival Time)")
print(f"Earliest Time = {shark_last_meal_time} + {SHARK_SURVIVAL_TIME} = {final_time} minutes.")
>>> 51