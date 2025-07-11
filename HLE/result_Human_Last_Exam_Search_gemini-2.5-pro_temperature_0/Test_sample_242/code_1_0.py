import random
import heapq

def run_truck_simulation(seed=42):
    """
    Runs a discrete-event simulation of the distribution center's truck unloading process.
    """
    random.seed(seed)

    # --- Simulation Parameters ---
    SIMULATION_DURATION = 8 * 60  # 480 minutes
    MEAN_ARRIVAL_INTERVAL = 7.0   # minutes
    MIN_SERVICE_TIME = 20.0       # minutes
    MAX_SERVICE_TIME = 35.0       # minutes
    TRAVEL_TIME_TO_DOCK = 1.0     # minute
    NUMBER_OF_DOCKS = 3

    # --- State Variables ---
    unloaded_truck_count = 0
    waiting_truck_queue = []
    event_heap = []  # Min-heap for (event_time, event_type, data)
    dock_free_times = [0.0] * NUMBER_OF_DOCKS

    # --- Helper Functions ---
    def schedule_arrival(current_time):
        """Schedules the next truck arrival event."""
        time_to_next_arrival = random.expovariate(1.0 / MEAN_ARRIVAL_INTERVAL)
        arrival_time = current_time + time_to_next_arrival
        if arrival_time <= SIMULATION_DURATION:
            heapq.heappush(event_heap, (arrival_time, 'TRUCK_ARRIVAL', {}))

    def schedule_departure(service_start_time, dock_idx):
        """Schedules a service completion (departure) event."""
        service_time = random.uniform(MIN_SERVICE_TIME, MAX_SERVICE_TIME)
        departure_time = service_start_time + TRAVEL_TIME_TO_DOCK + service_time
        dock_free_times[dock_idx] = departure_time
        heapq.heappush(event_heap, (departure_time, 'SERVICE_COMPLETION', {'dock_index': dock_idx}))

    # --- Simulation Initialization ---
    schedule_arrival(0)

    # --- Main Event Loop ---
    while event_heap:
        # Peek at the next event without removing it
        current_time, event_type, event_data = event_heap[0]

        # If the next event is past the end time, stop the simulation
        if current_time > SIMULATION_DURATION:
            break

        # Process the event
        heapq.heappop(event_heap)
        
        if event_type == 'TRUCK_ARRIVAL':
            schedule_arrival(current_time)
            
            free_dock_index = -1
            for i in range(NUMBER_OF_DOCKS):
                if dock_free_times[i] <= current_time:
                    free_dock_index = i
                    break
            
            if free_dock_index != -1:
                schedule_departure(current_time, free_dock_index)
            else:
                waiting_truck_queue.append(current_time)

        elif event_type == 'SERVICE_COMPLETION':
            unloaded_truck_count += 1
            completed_dock_index = event_data['dock_index']
            
            if waiting_truck_queue:
                waiting_truck_queue.pop(0)
                schedule_departure(current_time, completed_dock_index)

    return unloaded_truck_count

# Execute the simulation
unloaded_trucks = run_truck_simulation(seed=42)