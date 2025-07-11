import simpy
import random

def run_simulation():
    """
    This function sets up and runs the entire logistics simulation.
    """
    # --- Simulation Parameters ---
    RANDOM_SEED = 42          # For reproducibility
    SIM_TIME = 8 * 60         # Simulation time in minutes (8 hours)
    NUM_DOCKS = 3             # Number of docks
    MEAN_ARRIVAL = 7          # Mean time between truck arrivals (minutes)
    SERVICE_MIN = 20          # Minimum service time (minutes)
    SERVICE_MAX = 35          # Maximum service time (minutes)
    TRAVEL_TIME = 1           # Time for a truck to travel to the dock (minutes)

    # --- Statistics Counter ---
    # We use a list to hold the count so it can be modified by the truck function.
    unloaded_trucks_counter = [0]

    def truck(env, docks, counter):
        """
        Represents a truck's journey: arrival, queuing for a dock, unloading, and departure.
        """
        # 1. Truck requests a dock.
        # The 'with' statement handles requesting and releasing the resource automatically.
        with docks.request() as request:
            yield request  # Wait in the yard until a dock is available

            # 2. Travel to the now-available dock.
            yield env.timeout(TRAVEL_TIME)

            # 3. Perform the unloading service.
            service_duration = random.uniform(SERVICE_MIN, SERVICE_MAX)
            yield env.timeout(service_duration)

        # 4. Service is complete, truck has left the dock.
        # Increment the counter for unloaded trucks.
        counter[0] += 1

    def source(env, docks, counter):
        """Generates trucks arriving at the DC according to the exponential distribution."""
        while True:
            # Calculate the time until the next truck arrives.
            time_to_next_arrival = random.expovariate(1.0 / MEAN_ARRIVAL)
            yield env.timeout(time_to_next_arrival)
            
            # A new truck arrives, start its process.
            env.process(truck(env, docks, counter))

    # --- Main Simulation Setup ---
    # Set the random seed to ensure the result is the same every time.
    random.seed(RANDOM_SEED)

    # Create a simpy simulation environment.
    env = simpy.Environment()

    # Create the dock resource with a capacity of NUM_DOCKS.
    docks = simpy.Resource(env, capacity=NUM_DOCKS)

    # Start the truck arrival process.
    env.process(source(env, docks, unloaded_trucks_counter))

    # Run the simulation for the specified time.
    env.run(until=SIM_TIME)

    # --- Output the Final Result ---
    # The final equation shows the number of trucks unloaded in the given time.
    total_unloaded = unloaded_trucks_counter[0]
    hours_simulated = SIM_TIME / 60
    
    print(f"Total trucks unloaded in {hours_simulated:.0f} hours = {total_unloaded}")
    
# Execute the simulation
if __name__ == "__main__":
    run_simulation()