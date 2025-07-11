import random
import heapq

def run_dc_simulation():
    """
    Simulates truck arrivals and unloading at a Distribution Center for one 8-hour day.
    """
    # --- 1. Simulation Parameters ---
    SIMULATION_TIME_MINUTES = 8 * 60
    MEAN_ARRIVAL_MINUTES = 7.0
    NUM_DOCKS = 3
    MIN_SERVICE_MINUTES = 20.0
    MAX_SERVICE_MINUTES = 35.0
    TRAVEL_TIME_MINUTES = 1.0

    # Using a fixed seed for reproducible results
    random.seed(0)

    # --- 2. Generate Truck Arrivals ---
    arrival_times = []
    current_arrival_time = 0
    while True:
        # Inter-arrival time follows an exponential distribution
        inter_arrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_MINUTES)
        current_arrival_time += inter_arrival_time
        if current_arrival_time > SIMULATION_TIME_MINUTES:
            break
        arrival_times.append(current_arrival_time)

    # --- 3. Initialize Docks and Counters ---
    # Use a min-heap to keep track of when each dock becomes free.
    # Initially, all docks are free at time 0.
    dock_free_times = [0.0] * NUM_DOCKS
    heapq.heapify(dock_free_times)
    
    unloaded_trucks = 0

    # --- 4. Process Each Truck Arrival ---
    for arrival_time in arrival_times:
        # Get the earliest time a dock becomes available
        earliest_dock_free_time = heapq.heappop(dock_free_times)

        # A truck starts moving to the dock after it arrives AND the dock is free.
        start_move_time = max(arrival_time, earliest_dock_free_time)

        # Service begins after the travel time to the dock.
        service_start_time = start_move_time + TRAVEL_TIME_MINUTES

        # Service duration is uniform between 20 and 35 minutes.
        service_duration = random.uniform(MIN_SERVICE_MINUTES, MAX_SERVICE_MINUTES)

        # Calculate the final completion time.
        completion_time = service_start_time + service_duration

        # --- 5. Count if Unloaded Within the Day ---
        # The truck is counted only if its service finishes within the 8-hour shift.
        if completion_time <= SIMULATION_TIME_MINUTES:
            unloaded_trucks += 1
        
        # This dock is now occupied until the completion time. Add it back to the heap.
        heapq.heappush(dock_free_times, completion_time)

    # --- 6. Output the Results ---
    print("--- Distribution Center Simulation Results ---")
    print(f"Simulation Time: {SIMULATION_TIME_MINUTES} minutes")
    print(f"Number of Docks: {NUM_DOCKS}")
    print(f"Mean Truck Arrival: {MEAN_ARRIVAL_MINUTES} minutes (exponential)")
    print(f"Service Time: [{MIN_SERVICE_MINUTES}, {MAX_SERVICE_MINUTES}] minutes (uniform)")
    print(f"Yard-to-Dock Travel Time: {TRAVEL_TIME_MINUTES} minute")
    print("--------------------------------------------")
    print(f"Total trucks that arrived: {len(arrival_times)}")
    print(f"Total trucks unloaded within 8 hours: {unloaded_trucks}")

if __name__ == '__main__':
    run_dc_simulation()