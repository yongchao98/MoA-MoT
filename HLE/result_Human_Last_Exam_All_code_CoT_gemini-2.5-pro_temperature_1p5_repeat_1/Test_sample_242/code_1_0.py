import random

def simulate_dc_operations():
    """
    Simulates truck arrivals and unloading at a Distribution Center for an 8-hour day.
    """
    # --- 1. Simulation Parameters ---
    SIMULATION_TIME_MINUTES = 8 * 60
    MEAN_ARRIVAL_MINUTES = 7.0
    NUM_DOCKS = 3
    MIN_SERVICE_MINUTES = 20.0
    MAX_SERVICE_MINUTES = 35.0
    TRAVEL_TIME_MINUTES = 1.0

    # Use a fixed seed for reproducible results
    random.seed(42)

    # --- 2. Generate Truck Arrivals ---
    arrival_times = []
    current_time = 0.0
    while True:
        # Time until the next truck arrives, following an exponential distribution
        inter_arrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_MINUTES)
        current_time += inter_arrival_time
        
        # Stop generating arrivals after the 8-hour mark
        if current_time > SIMULATION_TIME_MINUTES:
            break
        arrival_times.append(current_time)

    # --- Simulation State ---
    # Stores the time when each dock will become free
    dock_free_time = [0.0] * NUM_DOCKS
    unloaded_trucks_count = 0

    # --- 3. Process Trucks Chronologically ---
    for truck_arrival_time in arrival_times:
        # Find the dock that becomes free the earliest
        earliest_free_dock_time = min(dock_free_time)
        earliest_free_dock_index = dock_free_time.index(earliest_free_dock_time)

        # A truck starts moving to the dock when it arrives or when the dock is free, whichever is later
        start_move_time = max(truck_arrival_time, earliest_free_dock_time)
        
        # Service can only begin after the truck travels from the yard to the dock
        service_start_time = start_move_time + TRAVEL_TIME_MINUTES
        
        # Generate a random service time from a uniform distribution
        service_duration = random.uniform(MIN_SERVICE_MINUTES, MAX_SERVICE_MINUTES)
        
        # Calculate when the service will be completed
        service_finish_time = service_start_time + service_duration
        
        # The dock is now occupied until the service finish time
        dock_free_time[earliest_free_dock_index] = service_finish_time
        
        # --- 4. Count Unloaded Trucks ---
        # A truck is counted as "unloaded" if its service finishes within the 8-hour day
        if service_finish_time <= SIMULATION_TIME_MINUTES:
            unloaded_trucks_count += 1
            
    # --- 5. Output the Result ---
    print(f"Total trucks arrived during the 8-hour day: {len(arrival_times)}")
    print(f"Total trucks unloaded within the 8-hour day: {unloaded_trucks_count}")

# Run the simulation
simulate_dc_operations()