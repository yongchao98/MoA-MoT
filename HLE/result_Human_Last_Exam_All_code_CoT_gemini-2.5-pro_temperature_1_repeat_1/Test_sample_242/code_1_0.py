import heapq
import random
from collections import deque

def simulate_dc_operations():
    """
    Simulates truck arrivals and service at a Distribution Center for one day.
    """
    # --- Simulation Parameters ---
    # Set a seed for reproducibility. The result will be the same on every run.
    random.seed(42)

    SIMULATION_HOURS = 8
    SIMULATION_MINUTES = SIMULATION_HOURS * 60  # 480 minutes
    NUM_DOCKS = 3
    MEAN_ARRIVAL_MINUTES = 7.0  # Exponential distribution
    MIN_SERVICE_MINUTES = 20    # Uniform distribution
    MAX_SERVICE_MINUTES = 35    # Uniform distribution
    TRAVEL_TIME_MINUTES = 1     # From yard to dock

    # --- Simulation State ---
    event_queue = []  # A min-heap to store future events: (time, event_type, truck_id)
    waiting_yard = deque()  # A FIFO queue for trucks waiting for a dock
    available_docks = NUM_DOCKS
    unloaded_truck_count = 0
    next_truck_id = 0

    # --- Event Types (as constants for clarity) ---
    TRUCK_ARRIVAL = "TRUCK_ARRIVAL"
    DOCK_ARRIVAL = "DOCK_ARRIVAL"
    SERVICE_COMPLETED = "SERVICE_COMPLETED"

    # --- Helper function to schedule an event ---
    def schedule_event(time, event_type, truck_id):
        heapq.heappush(event_queue, (time, event_type, truck_id))

    # --- Initialize Simulation ---
    # Schedule the very first truck arrival
    first_arrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_MINUTES)
    schedule_event(first_arrival_time, TRUCK_ARRIVAL, next_truck_id)
    next_truck_id += 1

    # --- Main Simulation Loop ---
    # The loop continues as long as there are events in the queue.
    while event_queue:
        current_time, event_type, truck_id = heapq.heappop(event_queue)

        # --- Event Handling Logic ---
        if event_type == TRUCK_ARRIVAL:
            # A new truck has arrived at the DC entrance.
            
            # If the arrival is within the working day, schedule the *next* arrival.
            if current_time <= SIMULATION_MINUTES:
                next_arrival_time = current_time + random.expovariate(1.0 / MEAN_ARRIVAL_MINUTES)
                # We only schedule future arrivals if they happen within the 8-hour window.
                if next_arrival_time <= SIMULATION_MINUTES:
                    schedule_event(next_arrival_time, TRUCK_ARRIVAL, next_truck_id)
                    next_truck_id += 1

            # Check for an available dock for the current truck.
            if available_docks > 0:
                available_docks -= 1
                schedule_event(current_time + TRAVEL_TIME_MINUTES, DOCK_ARRIVAL, truck_id)
            else:
                waiting_yard.append(truck_id)

        elif event_type == DOCK_ARRIVAL:
            # The truck has reached the dock, service begins.
            service_time = random.uniform(MIN_SERVICE_MINUTES, MAX_SERVICE_MINUTES)
            completion_time = current_time + service_time
            schedule_event(completion_time, SERVICE_COMPLETED, truck_id)

        elif event_type == SERVICE_COMPLETED:
            # The truck has finished unloading.
            # Only count if unloaded within the 8-hour day.
            if current_time <= SIMULATION_MINUTES:
                unloaded_truck_count += 1
            
            available_docks += 1

            # Check if any trucks are waiting in the yard to be serviced.
            if waiting_yard:
                next_in_line_id = waiting_yard.popleft()
                available_docks -= 1
                schedule_event(current_time + TRAVEL_TIME_MINUTES, DOCK_ARRIVAL, next_in_line_id)

    # --- Final Output ---
    print("Logistics DC Simulation")
    print("-" * 35)
    print("Parameters:")
    print(f"  Simulation Duration: {SIMULATION_HOURS} hours ({SIMULATION_MINUTES} minutes)")
    print(f"  Number of Docks: {NUM_DOCKS}")
    print(f"  Truck Arrival Mean: {MEAN_ARRIVAL_MINUTES} minutes (Exponential)")
    print(f"  Service Time: [{MIN_SERVICE_MINUTES}, {MAX_SERVICE_MINUTES}] minutes (Uniform)")
    print(f"  Yard to Dock Travel Time: {TRAVEL_TIME_MINUTES} minute")
    print("-" * 35)
    print("Result:")
    # The final "equation" is the final count based on the inputs.
    print(f"Total trucks unloaded in {SIMULATION_HOURS} hours = {unloaded_truck_count}")

if __name__ == '__main__':
    simulate_dc_operations()