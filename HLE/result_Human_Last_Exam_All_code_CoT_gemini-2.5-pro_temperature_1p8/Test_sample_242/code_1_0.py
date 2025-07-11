import random

def simulate_dc_operations():
    """
    Simulates truck arrivals and unloadings at a Distribution Center for an 8-hour day.
    """
    # --- Simulation Constants ---
    SIMULATION_DURATION = 8 * 60  # minutes
    NUM_DOCKS = 3
    MEAN_INTERARRIVAL_TIME = 7      # minutes, for exponential distribution
    SERVICE_TIME_MIN = 20           # minutes, for uniform distribution
    SERVICE_TIME_MAX = 35           # minutes, for uniform distribution
    TRUCK_TRAVEL_TIME = 1           # minute

    # Use a fixed seed for reproducible results
    random.seed(42)

    # --- Simulation State Variables ---
    # List to store the time when each dock becomes available
    docks_availability_times = [0.0] * NUM_DOCKS
    
    # Counter for trucks that finish unloading within the simulation duration
    trucks_unloaded_count = 0
    
    # Track the time of the current truck's arrival
    current_arrival_time = 0.0

    # --- Simulation Loop ---
    # The loop generates new truck arrivals until an arrival occurs after the 8-hour mark.
    while True:
        # Calculate the time for the next arrival
        interarrival_time = random.expovariate(1.0 / MEAN_INTERARRIVAL_TIME)
        current_arrival_time += interarrival_time

        # If the truck arrives after the simulation ends, stop the simulation
        if current_arrival_time > SIMULATION_DURATION:
            break

        # A truck has arrived. Find the earliest available dock.
        # The truck must go to the dock that has been waiting the longest (or is free).
        earliest_free_time = min(docks_availability_times)
        
        # Get the index of that dock to assign the truck to it.
        dock_index_for_truck = docks_availability_times.index(earliest_free_time)

        # The service can only start after the truck has arrived AND the dock is free.
        # Then, we add the 1-minute travel time from the yard to the dock.
        service_start_time = max(current_arrival_time, earliest_free_time) + TRUCK_TRAVEL_TIME
        
        # Determine the service duration for this truck.
        service_duration = random.uniform(SERVICE_TIME_MIN, SERVICE_TIME_MAX)
        
        # Calculate when the service will be completed.
        service_completion_time = service_start_time + service_duration
        
        # The assigned dock is now occupied and will become available at the completion time.
        docks_availability_times[dock_index_for_truck] = service_completion_time
        
        # We only count the truck if its unloading process is finished within the 8-hour day.
        if service_completion_time <= SIMULATION_DURATION:
            trucks_unloaded_count += 1
            
    print(f"Number of trucks unloaded in the 8-hour day: {trucks_unloaded_count}")

simulate_dc_operations()
<<<52>>>