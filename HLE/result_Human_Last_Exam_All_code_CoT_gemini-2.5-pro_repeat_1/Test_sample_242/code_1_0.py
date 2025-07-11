import simpy
import random

# --- Simulation Parameters ---
SIM_TIME = 8 * 60         # Simulation time in minutes (8 hours)
NUM_DOCKS = 3             # Number of docks for unloading
MEAN_ARRIVAL = 7          # Mean time between truck arrivals in minutes
MIN_SERVICE = 20          # Minimum unloading time in minutes
MAX_SERVICE = 35          # Maximum unloading time in minutes
TRAVEL_TIME = 1           # Travel time from yard to dock in minutes
RANDOM_SEED = 42          # Seed for reproducibility

def truck(env, name, docks, unloaded_counter):
    """
    Represents a single truck's journey: arriving, waiting for a dock,
    traveling to it, getting unloaded, and leaving.
    """
    # 1. Truck arrives at the DC
    arrival_time = env.now
    # print(f'{arrival_time:.2f}: {name} arrives at the DC.')

    # 2. Request a dock
    with docks.request() as request:
        yield request

        # 3. Check if waiting was required and add travel time
        time_got_dock = env.now
        # print(f'{time_got_dock:.2f}: {name} is assigned a dock.')
        if time_got_dock > arrival_time:
            # This truck had to wait in the yard, so it incurs travel time
            # print(f'{time_got_dock:.2f}: {name} must travel from yard to dock.')
            yield env.timeout(TRAVEL_TIME)
            # print(f'{env.now:.2f}: {name} has reached the dock.')

        # 4. Perform the unloading service
        service_time = random.uniform(MIN_SERVICE, MAX_SERVICE)
        # print(f'{env.now:.2f}: {name} starts unloading (will take {service_time:.2f} mins).')
        yield env.timeout(service_time)

        # 5. Unloading is complete
        # print(f'{env.now:.2f}: {name} has finished unloading.')
        unloaded_counter['count'] += 1

def source(env, docks, unloaded_counter):
    """Generates trucks arriving at the DC."""
    truck_id = 0
    while True:
        # Generate the time until the next truck arrives
        interarrival_time = random.expovariate(1.0 / MEAN_ARRIVAL)
        yield env.timeout(interarrival_time)
        
        # Create a new truck process
        truck_id += 1
        env.process(truck(env, f'Truck {truck_id}', docks, unloaded_counter))

if __name__ == '__main__':
    # --- Setup ---
    # Set the seed for random number generation to ensure reproducible results
    random.seed(RANDOM_SEED)

    # Create a simpy environment
    env = simpy.Environment()

    # Create the docks as a shared resource
    docks = simpy.Resource(env, capacity=NUM_DOCKS)

    # Create a counter to track completed unloadings
    unloaded_truck_count = {'count': 0}

    # --- Run Simulation ---
    # Start the truck arrival process
    env.process(source(env, docks, unloaded_truck_count))
    
    # Run the simulation for the specified duration
    env.run(until=SIM_TIME)

    # --- Output Results ---
    final_count = unloaded_truck_count['count']
    print(f"Simulation of an 8-hour day ({SIM_TIME} minutes)")
    print(f"Number of docks: {NUM_DOCKS}")
    print(f"Average truck arrival rate: every {MEAN_ARRIVAL} minutes")
    print(f"Service time: uniform between {MIN_SERVICE} and {MAX_SERVICE} minutes")
    print("-" * 40)
    print(f"Total number of trucks unloaded: {final_count}")
    print("<<<" + str(final_count) + ">>>")