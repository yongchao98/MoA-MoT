import simpy
import random

def run_simulation():
    """
    Sets up and runs the DC simulation.
    """
    # --- Simulation Parameters ---
    RANDOM_SEED = 42
    SIM_TIME = 8 * 60  # Simulation time in minutes
    NUM_DOCKS = 3
    MEAN_ARRIVAL = 7  # Mean time between truck arrivals (minutes)
    MIN_SERVICE = 20  # Min service time (minutes)
    MAX_SERVICE = 35  # Max service time (minutes)
    TRAVEL_TIME = 1   # Time for a truck to travel to the dock (minutes)

    # A simple class to hold simulation statistics
    class Stats:
        def __init__(self):
            self.unloaded_count = 0

    def truck(env, dc_docks, stats):
        """
        Represents the lifecycle of a single truck: requesting a dock,
        traveling to it, getting unloaded, and leaving.
        """
        # A truck requests a dock. It will wait if none are available.
        with dc_docks.request() as request:
            yield request

            # It takes time to travel to the assigned dock.
            yield env.timeout(TRAVEL_TIME)

            # The unloading process takes a variable amount of time.
            service_time = random.uniform(MIN_SERVICE, MAX_SERVICE)
            yield env.timeout(service_time)

        # The 'with' statement automatically releases the dock.
        # Now we can count this truck as successfully unloaded.
        stats.unloaded_count += 1

    def setup(env, num_docks, mean_arrival, stats):
        """
        Generates new trucks arriving at the DC.
        """
        # Create the resource representing the docks
        dc_docks = simpy.Resource(env, capacity=num_docks)

        # Continuously generate trucks until the simulation ends
        while True:
            # Wait for the next truck to arrive
            yield env.timeout(random.expovariate(1.0 / mean_arrival))
            
            # Create a new truck process once it arrives
            env.process(truck(env, dc_docks, stats))

    # --- Simulation Execution ---
    random.seed(RANDOM_SEED)
    env = simpy.Environment()
    simulation_stats = Stats()

    # Start the process that generates trucks
    env.process(setup(env, NUM_DOCKS, MEAN_ARRIVAL, simulation_stats))

    # Run the simulation for the specified duration
    env.run(until=SIM_TIME)

    # --- Output Results ---
    print("--- Simulation of a Distribution Center ---")
    print(f"Duration of one work day: {SIM_TIME / 60} hours ({SIM_TIME} minutes)")
    print(f"Number of docks for service: {NUM_DOCKS}")
    print(f"Truck arrival rate (exponential): 1 truck every {MEAN_ARRIVAL} minutes on average")
    print(f"Truck service time (uniform): between {MIN_SERVICE} and {MAX_SERVICE} minutes")
    print(f"Truck travel time from yard to dock: {TRAVEL_TIME} minute")
    print("\n--- Result ---")
    print(f"Total number of trucks unloaded in the 8-hour day: {simulation_stats.unloaded_count}")

# Run the main simulation function
run_simulation()