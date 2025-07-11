import simpy
import random

def run_simulation():
    """
    This function sets up and runs the entire logistics simulation.
    """
    # --- Simulation Parameters ---
    SIM_TIME = 8 * 60       # Simulation duration in minutes (8 hours)
    MEAN_ARRIVAL = 7.0      # Average time between truck arrivals (minutes)
    NUM_DOCKS = 3           # Number of service docks
    MIN_SERVICE = 20.0      # Minimum service time (minutes)
    MAX_SERVICE = 35.0      # Maximum service time (minutes)
    TRAVEL_TIME = 1.0       # Travel time from yard to dock (minutes)
    RANDOM_SEED = 42        # Use a seed for reproducible results

    # A counter to track completed trucks. Using a list to make it mutable
    # across different function scopes without using a global variable.
    unloaded_truck_counter = [0]

    def truck_process(env, docks, unloaded_counter):
        """
        Models the lifecycle of a single truck: arrival, queuing, service, departure.
        """
        # A dock is requested upon arrival. If all are busy, the truck waits.
        with docks.request() as request:
            yield request

            # Once a dock is granted, simulate the travel time to it.
            yield env.timeout(TRAVEL_TIME)

            # At the dock, simulate the unloading service time.
            service_time = random.uniform(MIN_SERVICE, MAX_SERVICE)
            yield env.timeout(service_time)

            # Unloading is complete, increment the counter.
            unloaded_counter[0] += 1

    def truck_generator(env, docks, unloaded_counter):
        """Generates trucks arriving at the DC according to the exponential distribution."""
        while True:
            # Wait for the next truck to arrive.
            yield env.timeout(random.expovariate(1.0 / MEAN_ARRIVAL))
            # Start the process for the newly arrived truck.
            env.process(truck_process(env, docks, unloaded_counter))

    # --- Simulation Execution ---
    random.seed(RANDOM_SEED)
    env = simpy.Environment()
    docks = simpy.Resource(env, capacity=NUM_DOCKS)
    env.process(truck_generator(env, docks, unloaded_truck_counter))
    env.run(until=SIM_TIME)

    # --- Final Output ---
    # The final "equation" showing all input numbers and the final result.
    print(
        f"Result: {unloaded_truck_counter[0]} trucks unloaded "
        f"(Inputs: {NUM_DOCKS} docks, {MEAN_ARRIVAL} min exponential arrival, "
        f"[{MIN_SERVICE}, {MAX_SERVICE}] min uniform service, "
        f"{TRAVEL_TIME} min travel, {SIM_TIME} min total simulation)"
    )

if __name__ == "__main__":
    run_simulation()