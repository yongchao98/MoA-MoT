import random
import heapq

def simulate_dc_operations():
    """
    Simulates truck arrivals and servicing at a Distribution Center for one 8-hour day.
    """
    # --- Simulation Parameters ---
    SIM_DURATION_MINUTES = 8 * 60
    NUM_DOCKS = 3
    MEAN_ARRIVAL_INTERVAL = 7.0
    MIN_SERVICE_TIME = 20.0
    MAX_SERVICE_TIME = 35.0
    TRUCK_TRAVEL_TIME = 1.0

    # Using a fixed seed for reproducible results
    RANDOM_SEED = 42
    random.seed(RANDOM_SEED)

    # --- Simulation State Variables ---
    # A priority queue for events, storing tuples of (time, event_type)
    event_queue = []
    # A list representing the time each of the 3 docks will become free
    docks_free_at_time = [0.0] * NUM_DOCKS
    # A FIFO queue for trucks waiting for a dock
    waiting_queue = []
    # The final count of trucks that finish unloading within the simulation time
    unloaded_truck_count = 0

    # --- Helper function to schedule arrivals ---
    def schedule_next_arrival(current_time):
        """Schedules a new ARRIVAL event if it falls within the 8-hour day."""
        interarrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_INTERVAL)
        next_arrival_time = current_time + interarrival_time
        if next_arrival_time <= SIM_DURATION_MINUTES:
            heapq.heappush(event_queue, (next_arrival_time, 'ARRIVAL'))

    # --- Start the simulation by scheduling the first arrival ---
    schedule_next_arrival(0)

    # --- Main Event Loop ---
    # The loop continues as long as there are events to process.
    while event_queue:
        current_time, event_type = heapq.heappop(event_queue)

        # We only count events that complete within the 8-hour window.
        # Any event occurring after the cutoff time is ignored.
        if current_time > SIM_DURATION_MINUTES:
            continue

        # --- Handle ARRIVAL Event ---
        if event_type == 'ARRIVAL':
            # A new truck arrived. Schedule the next one immediately.
            schedule_next_arrival(current_time)

            # Sort docks by free time to find the one that is available earliest.
            docks_free_at_time.sort()
            
            # Check if the earliest-available dock is free now.
            if docks_free_at_time[0] <= current_time:
                # A dock is free. Assign the truck to it.
                # Service starts at current_time.
                service_duration = random.uniform(MIN_SERVICE_TIME, MAX_SERVICE_TIME)
                total_dock_occupation_time = TRUCK_TRAVEL_TIME + service_duration
                completion_time = current_time + total_dock_occupation_time
                
                # Update this dock's availability time for the next truck.
                docks_free_at_time[0] = completion_time
                # Schedule the DEPARTURE event for when unloading is finished.
                heapq.heappush(event_queue, (completion_time, 'DEPARTURE'))
            else:
                # All docks are busy. The truck must wait in the yard.
                waiting_queue.append(current_time)

        # --- Handle DEPARTURE Event ---
        elif event_type == 'DEPARTURE':
            # A truck has finished unloading. Increment the counter.
            unloaded_truck_count += 1
            
            # The dock is now free. Check if any trucks are waiting.
            if waiting_queue:
                # A truck is waiting. Take the first one from the queue.
                waiting_queue.pop(0)
                
                # Service for the waiting truck starts now (at current_time).
                service_duration = random.uniform(MIN_SERVICE_TIME, MAX_SERVICE_TIME)
                total_dock_occupation_time = TRUCK_TRAVEL_TIME + service_duration
                completion_time = current_time + total_dock_occupation_time

                # Assign the newly-freed dock to this truck. We sort to ensure we update
                # the correct dock (the one whose free time was earliest).
                docks_free_at_time.sort()
                docks_free_at_time[0] = completion_time
                # Schedule this truck's DEPARTURE event.
                heapq.heappush(event_queue, (completion_time, 'DEPARTURE'))
                
    # The final result of the simulation
    print(unloaded_truck_count)

if __name__ == '__main__':
    simulate_dc_operations()