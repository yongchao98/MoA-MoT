import random
import heapq
from collections import deque

def run_logistics_simulation():
    """
    Simulates a truck unloading process at a Distribution Center (DC)
    to determine how many trucks are unloaded in an 8-hour day.
    """
    # --- Simulation Parameters ---
    # The operational management estimates are used here.
    SIMULATION_DURATION_MINUTES = 8 * 60  # An 8-hour day
    NUMBER_OF_DOCKS = 3                  # Three docks for service
    MEAN_ARRIVAL_INTERVAL_MINUTES = 7      # 7-minute exponential arrival
    SERVICE_TIME_MIN_MINUTES = 20        # Uniform service time minimum [20, 35]
    SERVICE_TIME_MAX_MINUTES = 35        # Uniform service time maximum [20, 35]
    TRUCK_TRAVEL_TIME_MINUTES = 1        # Travel time from entrance to dock area

    # Use a fixed seed for reproducible simulation results
    RANDOM_SEED = 42
    random.seed(RANDOM_SEED)

    # --- State Variables ---
    # The event_queue stores tuples: (event_time, event_type, truck_id).
    # It is a min-heap, so the next chronological event is always at the top.
    event_queue = []
    
    # A deque is an efficient list-like container for adding and removing from both ends.
    # It will store the IDs of trucks waiting for a dock.
    waiting_yard_queue = deque()
    
    # Tracks how many of the 3 docks are currently serving a truck.
    busy_docks_count = 0
    
    # The final metric: a counter for trucks that successfully finished unloading.
    unloaded_trucks_count = 0
    
    # Each truck gets a unique ID to be tracked throughout the simulation.
    truck_id_counter = 0

    # --- Helper function to schedule the next truck arrival ---
    def schedule_next_arrival(current_time):
        nonlocal truck_id_counter
        # Time to the next arrival follows an exponential distribution.
        # The rate (lambda) is 1/mean.
        time_to_next_arrival = random.expovariate(1.0 / MEAN_ARRIVAL_INTERVAL_MINUTES)
        arrival_time = current_time + time_to_next_arrival
        
        # We only schedule new arrivals that happen within the 8-hour workday.
        if arrival_time < SIMULATION_DURATION_MINUTES:
            # heapq.heappush adds an item to the heap, maintaining the heap property.
            heapq.heappush(event_queue, (arrival_time, "TRUCK_ENTRANCE_ARRIVAL", truck_id_counter))
            truck_id_counter += 1

    # --- Start the Simulation ---
    # Schedule the very first truck to start the process.
    schedule_next_arrival(0)
    
    # --- Main Simulation Loop ---
    # The loop continues as long as there are events scheduled to occur.
    while event_queue:
        # Get the next event chronologically using heapq.heappop.
        current_time, event_type, truck_id = heapq.heappop(event_queue)
        
        # Stop the simulation if the event is past the 8-hour mark.
        # Any service completions scheduled after this time will not be counted.
        if current_time >= SIMULATION_DURATION_MINUTES:
            break
            
        # --- Event Handling Logic ---
        
        if event_type == "TRUCK_ENTRANCE_ARRIVAL":
            # A new truck has arrived at the DC.
            # 1. Schedule the *next* truck's arrival to keep the simulation going.
            schedule_next_arrival(current_time)
            
            # 2. This truck now travels for 1 minute to the dock/yard area.
            #    Schedule its arrival there as a new event.
            time_at_dock_area = current_time + TRUCK_TRAVEL_TIME_MINUTES
            heapq.heappush(event_queue, (time_at_dock_area, "DOCK_AREA_ARRIVAL", truck_id))
            
        elif event_type == "DOCK_AREA_ARRIVAL":
            # The truck is now ready for service. Check for a free dock.
            if busy_docks_count < NUMBER_OF_DOCKS:
                # A dock is free. The truck begins service immediately.
                busy_docks_count += 1
                
                # Determine service time from a uniform distribution of [20, 35] minutes.
                service_time = random.uniform(SERVICE_TIME_MIN_MINUTES, SERVICE_TIME_MAX_MINUTES)
                completion_time = current_time + service_time
                
                # Schedule the 'SERVICE_COMPLETION' event for when this truck is done.
                heapq.heappush(event_queue, (completion_time, "SERVICE_COMPLETION", truck_id))
            else:
                # All docks are busy. The truck must wait in the yard.
                waiting_yard_queue.append(truck_id)
                
        elif event_type == "SERVICE_COMPLETION":
            # A truck has finished being unloaded.
            unloaded_trucks_count += 1
            
            # The dock is now free. Check if any truck is waiting in the yard.
            if waiting_yard_queue:
                # A truck is waiting. It starts service immediately at the freed dock.
                waiting_truck_id = waiting_yard_queue.popleft()
                
                # Determine the service time for this waiting truck.
                service_time = random.uniform(SERVICE_TIME_MIN_MINUTES, SERVICE_TIME_MAX_MINUTES)
                completion_time = current_time + service_time
                
                # Schedule the completion event for this new service.
                # Note: busy_docks_count does not change because one truck left and another took its place.
                heapq.heappush(event_queue, (completion_time, "SERVICE_COMPLETION", waiting_truck_id))
            else:
                # No trucks are waiting. A dock simply becomes free.
                busy_docks_count -= 1
                
    # --- Final Output ---
    # The prompt asks for the final number, which is our result.
    print(f"Simulation considering an 8-hour day:")
    print(f"- Mean truck arrival interval: {MEAN_ARRIVAL_INTERVAL_MINUTES} minutes")
    print(f"- Number of docks: {NUMBER_OF_DOCKS}")
    print(f"- Service time: [{SERVICE_TIME_MIN_MINUTES}, {SERVICE_TIME_MAX_MINUTES}] minutes")
    print("-" * 30)
    print(f"Total number of trucks that will be unloaded: {unloaded_trucks_count}")

# Execute the simulation function
run_logistics_simulation()
<<<48>>>