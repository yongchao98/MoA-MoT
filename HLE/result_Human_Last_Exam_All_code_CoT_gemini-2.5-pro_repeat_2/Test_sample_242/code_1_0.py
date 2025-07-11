import simpy
import random

# --- Simulation Parameters ---
# These are the "numbers" that define the simulation model.
SIM_TIME_HOURS = 8
MEAN_ARRIVAL_MINS = 7
NUM_DOCKS = 3
MIN_SERVICE_MINS = 20
MAX_SERVICE_MINS = 35
TRAVEL_TIME_MINS = 1
RANDOM_SEED = 42 # Use a seed for reproducible results

def truck_process(env, truck_id, docks, unloaded_truck_counter):
    """
    Represents the lifecycle of a single truck from arrival to departure.
    """
    # 1. Truck requests a dock
    with docks.request() as request:
        yield request  # Wait in the yard (queue) for a dock

        # 2. Travel to the dock
        yield env.timeout(TRAVEL_TIME_MINS)

        # 3. Get unloaded at the dock
        service_time = random.uniform(MIN_SERVICE_MINS, MAX_SERVICE_MINS)
        yield env.timeout(service_time)

    # 4. Unloading is complete, increment the counter
    unloaded_truck_counter[0] += 1

def truck_generator(env, docks, unloaded_truck_counter):
    """
    Generates new trucks at intervals following an exponential distribution.
    """
    truck_count = 0
    while True:
        # Wait for the next truck to arrive
        interarrival_time = random.expovariate(1.0 / MEAN_ARRIVAL_MINS)
        yield env.timeout(interarrival_time)

        # A new truck has arrived, start its process
        truck_count += 1
        env.process(truck_process(env, f"Truck {truck_count}", docks, unloaded_truck_counter))

def run_simulation():
    """
    Sets up and runs the simulation, then prints the results.
    You may need to install simpy: pip install simpy
    """
    # Set the random seed
    random.seed(RANDOM_SEED)

    # Use a list to hold the counter, making it mutable across functions
    unloaded_trucks_count = [0]

    # Create a simpy environment
    env = simpy.Environment()

    # Create the dock resources
    docks = simpy.Resource(env, capacity=NUM_DOCKS)

    # Start the truck generator process
    env.process(truck_generator(env, docks, unloaded_trucks_count))

    # Run the simulation for the specified duration
    simulation_duration_mins = SIM_TIME_HOURS * 60
    env.run(until=simulation_duration_mins)

    # --- Output the results ---
    print("--- Simulation Parameters (The numbers in the 'equation') ---")
    print(f"Simulation duration: {simulation_duration_mins} minutes ({SIM_TIME_HOURS} hours)")
    print(f"Number of docks for service: {NUM_DOCKS}")
    print(f"Truck arrival rate (exponential mean): {MEAN_ARRIVAL_MINS} minutes")
    print(f"Truck service time (uniform range): [{MIN_SERVICE_MINS}, {MAX_SERVICE_MINS}] minutes")
    print(f"Truck travel time from yard to dock: {TRAVEL_TIME_MINS} minute")
    print("-" * 55)
    print("--- Simulation Result ---")
    print(f"Total number of trucks unloaded: {unloaded_trucks_count[0]}")


if __name__ == '__main__':
    run_simulation()