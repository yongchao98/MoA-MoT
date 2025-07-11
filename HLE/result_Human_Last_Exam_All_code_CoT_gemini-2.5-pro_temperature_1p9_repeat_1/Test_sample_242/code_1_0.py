import random

def simulate_dc_operations():
    """
    Simulates truck arrivals and unloading at a Distribution Center (DC)
    for an 8-hour day to determine the number of trucks unloaded.
    """
    # --- 1. Define Simulation Parameters ---
    # Using a fixed seed for the random number generator to ensure reproducibility.
    random.seed(42)

    SIMULATION_DURATION_MINS = 8 * 60
    MEAN_ARRIVAL_INTERVAL_MINS = 7
    NUM_DOCKS = 3
    SERVICE_TIME_MIN_MINS = 20
    SERVICE_TIME_MAX_MINS = 35
    TRAVEL_TIME_MINS = 1

    # --- 2. Generate Truck Arrivals ---
    # Generate a list of all truck arrival times within the 8-hour window.
    arrival_times = []
    current_time = 0.0
    while True:
        # Time to next arrival follows an exponential distribution.
        interarrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_INTERVAL_MINS)
        current_time += interarrival_time
        if current_time > SIMULATION_DURATION_MINS:
            break
        arrival_times.append(current_time)

    # --- 3. Simulate the Unloading Process ---
    # This list tracks the time when each dock becomes free. Initialize all to 0.
    docks_free_at = [0.0] * NUM_DOCKS
    completion_times = []

    # Process each truck in the order of its arrival.
    for truck_arrival_time in arrival_times:
        # Find the dock that will be free earliest.
        # This models the truck going to the first available dock.
        earliest_free_time = min(docks_free_at)
        earliest_free_dock_index = docks_free_at.index(earliest_free_time)

        # Service can only start after the truck arrives AND the dock is free.
        # This line implicitly handles the waiting queue.
        service_start_time = max(truck_arrival_time, earliest_free_time)

        # Calculate the random service time for this truck (uniform distribution).
        service_duration = random.uniform(SERVICE_TIME_MIN_MINS, SERVICE_TIME_MAX_MINS)

        # Calculate when the service for this truck will be completed.
        service_completion_time = service_start_time + TRAVEL_TIME_MINS + service_duration

        # Record this completion time to count later.
        completion_times.append(service_completion_time)

        # Update the availability of the chosen dock to the new completion time.
        docks_free_at[earliest_free_dock_index] = service_completion_time

    # --- 4. Count Trucks Unloaded Within the 8-Hour Day ---
    # A truck is counted if its unloading process finishes within the workday.
    unloaded_truck_count = 0
    for time in completion_times:
        if time <= SIMULATION_DURATION_MINS:
            unloaded_truck_count += 1
            
    # --- 5. Print the Final Result ---
    print(f"Simulation of a Distribution Center")
    print("-" * 40)
    print(f"Parameters:")
    print(f"  Simulation Duration: {SIMULATION_DURATION_MINS} minutes (8 hours)")
    print(f"  Number of Docks: {NUM_DOCKS}")
    print(f"  Mean Truck Arrival Interval: {MEAN_ARRIVAL_INTERVAL_MINS} minutes")
    print(f"  Service Time: Uniform({SERVICE_TIME_MIN_MINS}, {SERVICE_TIME_MAX_MINS}) minutes")
    print(f"  Yard-to-Dock Travel Time: {TRAVEL_TIME_MINS} minute")
    print("-" * 40)
    print(f"Results:")
    print(f"  Total trucks that arrived during the day: {len(arrival_times)}")
    print(f"  Number of trucks unloaded = {unloaded_truck_count}")


if __name__ == '__main__':
    simulate_dc_operations()