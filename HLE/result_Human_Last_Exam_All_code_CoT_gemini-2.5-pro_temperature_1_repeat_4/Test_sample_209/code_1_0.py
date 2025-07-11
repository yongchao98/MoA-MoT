# This script simulates the optimal strategy to find the earliest time the pond is empty.

# --- Initial State ---
pond_goldfish = 10
basket_goldfish = 0
shark_meal_times = []
time = 0

print("--- Phase 1: Fisherman uses Rod B until basket is full ---")
print(f"Time {time:2}: Start. Pond has {pond_goldfish} goldfish.")

# The fisherman uses Rod B, catching a fish every 5 minutes.
# The shark eats a fish every 10 minutes.
# We check events at 5-minute intervals up to t=20.
for i in range(1, 5):
    catch_time = i * 5
    time = catch_time
    
    # Fisherman's catch
    pond_goldfish -= 1
    basket_goldfish += 1
    print(f"Time {time:2}: Fisherman catches a goldfish. Pond: {pond_goldfish}, Basket: {basket_goldfish}.")

    # Shark's meal (occurs at t=10, 20, 30...)
    if time % 10 == 0:
        pond_goldfish -= 1
        shark_meal_times.append(time)
        print(f"Time {time:2}: Shark eats a goldfish. Pond: {pond_goldfish}. Shark has now eaten {len(shark_meal_times)} time(s).")

# --- Phase 2: Fisherman gets a new basket ---
print("\n--- Phase 2: Basket is full, Fisherman leaves for 20 minutes ---")
time_of_departure = 20
trip_duration = 20
time_of_return = time_of_departure + trip_duration
print(f"Time {time_of_departure:2}: Fisherman leaves. He will return at t={time_of_return}.")

# The shark continues to eat while the fisherman is away.
last_meal = shark_meal_times[-1]
next_meal_time = last_meal + 10
if next_meal_time < time_of_return:
    time = next_meal_time
    pond_goldfish -= 1
    shark_meal_times.append(time)
    print(f"Time {time:2}: Shark eats a goldfish. Pond: {pond_goldfish}. Shark has now eaten {len(shark_meal_times)} time(s).")

# The shark's next meal is at t=40, exactly when the fisherman returns.
last_meal = shark_meal_times[-1]
next_meal_time = last_meal + 10
if next_meal_time == time_of_return:
    time = next_meal_time
    pond_goldfish -= 1
    shark_meal_times.append(time)
    print(f"Time {time:2}: Shark eats a goldfish. Pond: {pond_goldfish}. Shark has now eaten {len(shark_meal_times)} time(s).")

time = time_of_return
basket_goldfish = 0 # He has a new, empty basket.
print(f"Time {time:2}: Fisherman returns. The pond has {pond_goldfish} goldfish.")


# --- Phase 3: Fisherman clears the remaining goldfish ---
print("\n--- Phase 3: Clearing the rest of the pond ---")
while pond_goldfish > 0:
    time += 5 # 5 minutes per fish with Rod B
    pond_goldfish -= 1
    basket_goldfish += 1
    print(f"Time {time:2}: Fisherman catches a goldfish. Pond: {pond_goldfish}, Basket: {basket_goldfish}.")

pond_empty_time = time
print(f"\nTime {pond_empty_time:2}: All goldfish are gone from the pond.")


# --- Phase 4: Calculate the shark's death ---
print("\n--- Phase 4: Determining when the shark dies ---")
last_shark_meal = shark_meal_times[-1]
print(f"The shark's last meal was at t={last_shark_meal}. It cannot eat at t={last_shark_meal + 10} because the pond is empty.")
print("The shark's 11-minute starvation timer starts from its last meal.")

# Check for survival time extension
meals_in_last_60_min = 0
for meal_t in shark_meal_times:
    if pond_empty_time - 60 < meal_t <= pond_empty_time:
        meals_in_last_60_min += 1

print(f"In the 60 minutes before the pond was cleared of goldfish, the shark ate {meals_in_last_60_min} times.")

starvation_duration = 11
if meals_in_last_60_min > 4:
    starvation_duration = 11 + 4
    print("Because it ate more than 4 fish, its survival time is extended by 4 minutes to 15 minutes.")
else:
    print("Because it did not eat more than 4 fish, its survival time is the standard 11 minutes.")

shark_death_time = last_shark_meal + starvation_duration

print("\nThe final time is when the shark dies. This is calculated from its last meal.")
print("\nFinal Calculation:")
print(f"Time of last meal + Starvation duration = Total Time")
print(f"{last_shark_meal} minutes + {starvation_duration} minutes = {shark_death_time} minutes")

print(f"\nThe earliest time the pond is completely empty of fish is {shark_death_time} minutes.")