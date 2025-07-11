import random
import heapq
import collections

def run_dc_simulation():
    """
    Simulates a Distribution Center (DC) for an 8-hour day to determine
    the number of trucks unloaded.
    """
    # 1. Define Simulation Parameters
    SIM_DURATION_MINUTES = 8 * 60
    NUM_DOCKS = 3
    MEAN_ARRIVAL_MINUTES = 7.0
    SERVICE_TIME_MIN_MINUTES = 20.0
    SERVICE_TIME_MAX_MINUTES = 35.0
    TRAVEL_TIME_MINUTES = 1.0
    
    # Set a seed for reproducibility of results
    random.seed(42)

    # 2. Initialize State and Event Queue
    event_queue = []  # A min-heap for (time, event_type)
    waiting_queue = collections.deque() # Queue for trucks waiting for a dock
    
    current_time = 0.0
    available_docks = NUM_DOCKS
    unloaded_trucks_count = 0

    def start_service(time_of_dock_assignment):
        """Schedules a service completion event for a truck."""
        nonlocal available_docks, event_queue
        
        available_docks -= 1
        
        # Calculate when service starts (after travel) and ends
        service_start_time = time_of_dock_assignment + TRAVEL_TIME_MINUTES
        service_duration = random.uniform(SERVICE_TIME_MIN_MINUTES, SERVICE_TIME_MAX_MINUTES)
        service_end_time = service_start_time + service_duration
        
        # Add the 'SERVICE_END' event to the queue
        heapq.heappush(event_queue, (service_end_time, 'SERVICE_END'))

    # 3. Kickstart the simulation by scheduling the first arrival
    first_arrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_MINUTES)
    heapq.heappush(event_queue, (first_arrival_time, 'ARRIVAL'))

    # 4. Run the main simulation loop
    while event_queue:
        # Get the next event from the queue
        event_time, event_type = heapq.heappop(event_queue)
        
        # If the event happens after the simulation ends, stop
        if event_time > SIM_DURATION_MINUTES:
            break
        
        # Advance the simulation clock
        current_time = event_time

        if event_type == 'ARRIVAL':
            # Schedule the next truck's arrival
            next_arrival_time = current_time + random.expovariate(1.0 / MEAN_ARRIVAL_MINUTES)
            heapq.heappush(event_queue, (next_arrival_time, 'ARRIVAL'))

            # Process the current truck
            if available_docks > 0:
                start_service(current_time)
            else:
                waiting_queue.append(1) # Add a truck to the waiting yard
        
        elif event_type == 'SERVICE_END':
            unloaded_trucks_count += 1
            available_docks += 1

            # Check if there's a truck waiting in the yard for this freed dock
            if waiting_queue:
                waiting_queue.popleft()
                start_service(current_time) # The truck starts its journey now
    
    # 5. Output the results
    print("--- DC Simulation Parameters ---")
    print(f"Simulation duration: {SIM_DURATION_MINUTES / 60} hours ({SIM_DURATION_MINUTES} minutes)")
    print(f"Truck arrival rate (mean): 1 every {MEAN_ARRIVAL_MINUTES} minutes (exponential)")
    print(f"Number of docks for service: {NUM_DOCKS}")
    print(f"Truck unloading time (service): [{int(SERVICE_TIME_MIN_MINUTES)}, {int(SERVICE_TIME_MAX_MINUTES)}] minutes (uniform)")
    print(f"Truck travel time to dock: {int(TRAVEL_TIME_MINUTES)} minute")
    print("\n--- Simulation Result ---")
    print(f"Total trucks unloaded in one 8-hour day: {unloaded_trucks_count}")

# Run the simulation and print the results
if __name__ == "__main__":
    run_dc_simulation()
    
<<<68>>>