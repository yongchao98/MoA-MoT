import random

def simulate_dc_operations():
    """
    Simulates truck arrivals and unloading at a Distribution Center for one day.
    """
    # --- Simulation Parameters ---
    SIMULATION_DURATION = 8 * 60  # in minutes
    NUM_DOCKS = 3
    MEAN_ARRIVAL_TIME = 7.0       # minutes, exponential distribution
    SERVICE_TIME_MIN = 20         # minutes, uniform distribution
    SERVICE_TIME_MAX = 35         # minutes, uniform distribution
    TRAVEL_TIME = 1               # minute

    # Set a seed for the random number generator for reproducible results
    random.seed(42)

    # --- Simulation State Variables ---
    # A list to store the time when each dock becomes free.
    # Initially, all docks are free at time 0.
    docks_free_at_time = [0] * NUM_DOCKS

    unloaded_trucks_count = 0
    current_arrival_time = 0

    # --- Main Simulation Loop ---
    # The loop continues as long as new trucks arrive within the workday.
    while True:
        # Calculate the time until the next truck arrives using an exponential distribution.
        inter_arrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_TIME)
        
        # Advance the simulation time to this new truck's arrival.
        current_arrival_time += inter_arrival_time

        # If the truck arrives after the 8-hour day is over, no more trucks can be processed.
        if current_arrival_time > SIMULATION_DURATION:
            break

        # --- Process the Arrived Truck ---
        
        # Find which dock will be free at the earliest time.
        earliest_dock_free_time = min(docks_free_at_time)
        
        # Get the index of that dock.
        available_dock_index = docks_free_at_time.index(earliest_dock_free_time)

        # A truck can only start its journey to the dock after it arrives AND the dock is free.
        # So, the start time is the later of these two events.
        start_process_time = max(current_arrival_time, earliest_dock_free_time)

        # Generate a random service time for this specific truck from a uniform distribution.
        service_duration = random.uniform(SERVICE_TIME_MIN, SERVICE_TIME_MAX)

        # Calculate the final time the service will be completed.
        # This includes the time to start, the travel time, and the unloading time.
        service_finish_time = start_process_time + TRAVEL_TIME + service_duration

        # The chosen dock is now occupied and will become free at the service_finish_time.
        docks_free_at_time[available_dock_index] = service_finish_time

        # A truck is only counted if its unloading is finished within the 8-hour workday.
        if service_finish_time <= SIMULATION_DURATION:
            unloaded_trucks_count += 1
            
    # --- Output Results ---
    print("--- Simulation Parameters ---")
    print(f"Simulation duration: {SIMULATION_DURATION} minutes ({SIMULATION_DURATION/60} hours)")
    print(f"Number of docks: {NUM_DOCKS}")
    print(f"Mean truck arrival time: {MEAN_ARRIVAL_TIME} minutes")
    print(f"Service time range: [{SERVICE_TIME_MIN}, {SERVICE_TIME_MAX}] minutes")
    print(f"Travel time to dock: {TRAVEL_TIME} minute")
    print("\n--- Simulation Result ---")
    print(f"Total number of trucks unloaded = {unloaded_trucks_count}")

if __name__ == '__main__':
    simulate_dc_operations()
<<<53>>>