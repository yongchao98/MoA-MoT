import random
import heapq
from collections import deque

def run_logistics_simulation():
    """
    Simulates a truck unloading process at a Distribution Center for an 8-hour day.
    """
    # --- Simulation Parameters ---
    SIM_DURATION = 8 * 60      # Simulation time in minutes
    NUM_DOCKS = 3              # Number of available docks
    MEAN_ARRIVAL_TIME = 7.0    # Mean time between truck arrivals (minutes)
    SERVICE_TIME_MIN = 20.0    # Minimum unload time (minutes)
    SERVICE_TIME_MAX = 35.0    # Maximum unload time (minutes)
    TRAVEL_TIME = 1.0          # Time for a truck to travel to a dock (minute)
    RANDOM_SEED = 42           # Seed for reproducibility

    # Set the seed for random number generation
    random.seed(RANDOM_SEED)

    # --- State Variables ---
    # Event queue stores tuples: (event_time, event_type, data_dictionary)
    event_queue = []
    # Queue for trucks waiting in the yard
    waiting_trucks = deque()
    # List to track the time each dock becomes free
    docks_free_time = [0.0] * NUM_DOCKS
    
    unloaded_trucks_count = 0
    truck_id_counter = 0

    # --- Helper function to schedule events ---
    def schedule_arrival(current_time):
        nonlocal truck_id_counter
        # Calculate next arrival time using an exponential distribution
        interarrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_TIME)
        arrival_time = current_time + interarrival_time
        
        # Schedule the new arrival if it's within the simulation duration
        truck_id_counter += 1
        heapq.heappush(event_queue, (arrival_time, 'ARRIVAL', {'truck_id': truck_id_counter}))

    def schedule_unload(current_time, truck_id, dock_id):
        # Calculate when the unload will finish
        # Travel starts at current_time, unload starts after travel
        service_start_time = current_time + TRAVEL_TIME
        # Service time is uniformly distributed
        service_duration = random.uniform(SERVICE_TIME_MIN, SERVICE_TIME_MAX)
        service_end_time = service_start_time + service_duration
        
        # Mark the dock as busy until the service ends
        docks_free_time[dock_id] = service_end_time
        
        # Schedule the UNLOAD_END event
        heapq.heappush(event_queue, (service_end_time, 'UNLOAD_END', {'truck_id': truck_id, 'dock_id': dock_id}))

    # --- Start the simulation ---
    # Schedule the very first truck arrival at time 0
    schedule_arrival(0)
    # The first truck ID will be 1, but its arrival is after time 0.
    # We remove the initial dummy event used for scheduling.
    heapq.heappop(event_queue)
    truck_id_counter = 0
    schedule_arrival(0)

    # --- Main Simulation Loop ---
    # Process events as long as there are any and they are within the simulation time
    while event_queue and event_queue[0][0] <= SIM_DURATION:
        
        # Get the next event from the queue
        event_time, event_type, data = heapq.heappop(event_queue)
        
        # --- Event Handling Logic ---
        if event_type == 'ARRIVAL':
            truck_id = data['truck_id']
            
            # Schedule the *next* truck's arrival
            schedule_arrival(event_time)
            
            # Check for an available dock
            # Docks are sorted by their free time to always pick the one that is free first
            docks_free_time.sort()
            if docks_free_time[0] <= event_time:
                # A dock is free, send the truck to it
                schedule_unload(event_time, truck_id, dock_id=0)
            else:
                # All docks are busy, truck must wait in the yard
                waiting_trucks.append(truck_id)

        elif event_type == 'UNLOAD_END':
            dock_id = data['dock_id']
            
            # Increment count of unloaded trucks
            unloaded_trucks_count += 1
            
            # The dock is now free. Check if a truck is waiting in the yard.
            if waiting_trucks:
                # Get the next truck from the queue
                waiting_truck_id = waiting_trucks.popleft()
                # Send this truck to the now-free dock
                # The truck starts its travel from the yard at the current event_time
                schedule_unload(event_time, waiting_truck_id, dock_id)

    # --- Final Output ---
    print(f"Total trucks unloaded in an 8-hour day: {unloaded_trucks_count}")

if __name__ == '__main__':
    run_logistics_simulation()