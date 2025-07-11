import random

def run_dc_simulation():
    """
    Simulates truck arrivals and unloading at a distribution center for one day.
    """
    # --- 1. Simulation Parameters ---
    SIMULATION_HOURS = 8
    SIMULATION_MINUTES = SIMULATION_HOURS * 60
    NUM_DOCKS = 3
    MEAN_INTER_ARRIVAL_TIME = 7  # Exponential distribution mean
    MIN_SERVICE_TIME = 20          # Uniform distribution min
    MAX_SERVICE_TIME = 35          # Uniform distribution max
    TRAVEL_TIME = 1

    # To ensure the simulation gives the same result every time, we set a seed.
    random.seed(0)

    print("--- Simulation Parameters ---")
    print(f"Duration: {SIMULATION_MINUTES} minutes ({SIMULATION_HOURS} hours)")
    print(f"Number of Docks: {NUM_DOCKS}")
    print(f"Mean Truck Inter-arrival Time: {MEAN_INTER_ARRIVAL_TIME} minutes (Exponential)")
    print(f"Service Time: [{MIN_SERVICE_TIME}, {MAX_SERVICE_TIME}] minutes (Uniform)")
    print(f"Truck Travel Time to Dock: {TRAVEL_TIME} minute")
    print("---------------------------\n")

    # --- 2. Generate Truck Arrival Times ---
    arrival_times = []
    current_time = 0.0
    while True:
        # Calculate time until the next truck arrives
        inter_arrival_time = random.expovariate(1.0 / MEAN_INTER_ARRIVAL_TIME)
        current_time += inter_arrival_time

        if current_time > SIMULATION_MINUTES:
            break
        arrival_times.append(current_time)

    # --- 3. & 4. Simulate the Process and Count Unloaded Trucks ---
    # This list will hold the time at which each dock becomes free.
    dock_availability_times = [0.0] * NUM_DOCKS
    unloaded_truck_count = 0

    # Process each truck that arrives during the workday
    for arrival_time in arrival_times:
        # Find the earliest time a dock is free and which dock it is
        earliest_free_time = min(dock_availability_times)
        dock_index = dock_availability_times.index(earliest_free_time)

        # A truck can only start moving to the dock after it has arrived AND a dock is free.
        start_travel_time = max(arrival_time, earliest_free_time)
        
        # Service begins after the 1-minute travel to the dock.
        service_start_time = start_travel_time + TRAVEL_TIME

        # Determine the service duration for this specific truck.
        service_duration = random.uniform(MIN_SERVICE_TIME, MAX_SERVICE_TIME)

        # Calculate when the service for this truck will be completed.
        service_end_time = service_start_time + service_duration

        # This dock is now occupied until the service_end_time.
        # This must be updated regardless of whether the truck finishes on time,
        # as it affects the queue for subsequent trucks.
        dock_availability_times[dock_index] = service_end_time

        # Only count the truck if its unloading is completed within the 8-hour day.
        if service_end_time <= SIMULATION_MINUTES:
            unloaded_truck_count += 1

    # --- 5. Output the Result ---
    print(f"Total trucks arrived: {len(arrival_times)}")
    print(f"Result: A total of {unloaded_truck_count} trucks will be unloaded in an 8-hour day.")

if __name__ == "__main__":
    run_dc_simulation()

<<<48>>>