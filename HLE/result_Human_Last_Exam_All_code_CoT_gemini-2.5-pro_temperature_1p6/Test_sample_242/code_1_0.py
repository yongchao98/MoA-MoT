import simpy
import random

def run_simulation():
    """
    This function sets up and runs the logistics simulation.
    """
    # --- Simulation Parameters ---
    MEAN_ARRIVAL_TIME = 7      # Average minutes between truck arrivals (exponential)
    MIN_SERVICE_TIME = 20      # Minimum minutes for unloading (uniform)
    MAX_SERVICE_TIME = 35      # Maximum minutes for unloading (uniform)
    NUM_DOCKS = 3              # Number of available service docks
    TRAVEL_TIME = 1            # Minutes for a truck to travel from yard to dock
    SIMULATION_DURATION = 8 * 60 # Total simulation time in minutes (8-hour day)
    RANDOM_SEED = 42           # Seed for random number generator for reproducibility

    # --- Simulation State ---
    # We use a dictionary to hold the results, as it is a mutable object
    # that can be modified by the simulation processes.
    stats = {'unloaded_trucks': 0}

    def truck_process(env, name, docks, stats):
        """
        Represents the lifecycle of a single truck.
        It arrives, requests a dock, travels to it, gets unloaded, and leaves.
        """
        # A truck "arrives" and immediately requests a dock. Time spent here is queue time.
        with docks.request() as request:
            yield request  # Wait in the queue for a free dock

            # Once a dock is free, the truck starts traveling to it.
            yield env.timeout(TRAVEL_TIME)

            # The truck has reached the dock, and unloading begins.
            service_duration = random.uniform(MIN_SERVICE_TIME, MAX_SERVICE_TIME)
            yield env.timeout(service_duration)

            # Unloading is complete. Increment the counter.
            # This line is only reached if the service completes before the simulation ends.
            stats['unloaded_trucks'] += 1

    def truck_generator(env, docks, stats):
        """
        Generates new trucks arriving at the distribution center.
        """
        truck_id = 0
        while True:
            # Calculate the time until the next truck arrives.
            time_to_next_arrival = random.expovariate(1.0 / MEAN_ARRIVAL_TIME)
            yield env.timeout(time_to_next_arrival)

            # A new truck has arrived; create its process.
            truck_id += 1
            env.process(truck_process(env, f'Truck {truck_id}', docks, stats))

    # --- Simulation Execution ---
    random.seed(RANDOM_SEED)
    env = simpy.Environment()
    docks = simpy.Resource(env, capacity=NUM_DOCKS)
    env.process(truck_generator(env, docks, stats))
    env.run(until=SIMULATION_DURATION)

    # --- Output Results ---
    print(f"Simulation completed for an 8-hour day.")
    print(f"Number of docks: {NUM_DOCKS}")
    print(f"Arrival rate (mean): {MEAN_ARRIVAL_TIME} minutes")
    print(f"Service time (uniform): [{MIN_SERVICE_TIME}, {MAX_SERVICE_TIME}] minutes")
    print("---")
    # Final answer based on the simulation run.
    print(f"Total trucks unloaded in 8 hours: {stats['unloaded_trucks']}")

if __name__ == '__main__':
    run_simulation()