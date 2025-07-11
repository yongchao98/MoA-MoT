import simpy
import random

def solve_logistics_simulation():
    """
    Solves the logistics distribution center simulation problem.
    """

    # --- 1. Setup the Simulation Parameters ---
    SIM_TIME = 8 * 60      # Simulation time in minutes (8 hours)
    NUM_DOCKS = 3          # Number of docks for service
    MEAN_ARRIVAL = 7       # Average time between truck arrivals (minutes)
    MIN_SERVICE = 20       # Minimum service time (minutes)
    MAX_SERVICE = 35       # Maximum service time (minutes)
    TRAVEL_TIME = 1        # Travel time from yard to dock (minutes)
    RANDOM_SEED = 42       # Use a seed for reproducibility of results

    # This list is used as a mutable counter to track unloaded trucks
    unloaded_trucks_count = [0]

    # --- 2. Define the Simulation Processes ---

    def truck(env, name, docks, counter):
        """
        Represents the lifecycle of a single truck from arrival to departure.
        A truck arrives, requests a dock, travels to it, gets unloaded, and leaves.
        """
        # A truck requests one of the docks. This is where it waits in the yard if needed.
        with docks.request() as request:
            yield request # Wait in the queue until a dock is granted

            # Once a dock is assigned, the truck travels for 1 minute to reach it.
            yield env.timeout(TRAVEL_TIME)
            
            # The truck is now at the dock, and unloading begins.
            # We generate a random service time from a uniform distribution.
            service_time = random.uniform(MIN_SERVICE, MAX_SERVICE)
            yield env.timeout(service_time)
            
            # Unloading is complete, so we increment the counter.
            counter[0] += 1

    def source(env, docks, counter):
        """
        Generates new trucks arriving at the DC according to an exponential distribution.
        """
        truck_id = 0
        while True:
            # Generate the time until the next truck arrives.
            inter_arrival_time = random.expovariate(1.0 / MEAN_ARRIVAL)
            yield env.timeout(inter_arrival_time)
            
            truck_id += 1
            # Create a new truck process when a truck arrives.
            env.process(truck(env, f'Truck {truck_id}', docks, counter))

    # --- 3. Run the Simulation ---

    # Set the seed for the random number generator for consistent results.
    random.seed(RANDOM_SEED)

    # Create a simpy simulation environment.
    env = simpy.Environment()

    # Create the dock resources with a capacity of NUM_DOCKS.
    docks = simpy.Resource(env, capacity=NUM_DOCKS)

    # Start the truck arrival generator process.
    env.process(source(env, docks, unloaded_trucks_count))

    # Run the simulation for the specified duration (480 minutes).
    env.run(until=SIM_TIME)

    # --- 4. Report the Final Result ---
    
    # The final equation is the total number of trucks unloaded during the 8-hour day.
    # The code outputs each number involved in the final result statement.
    print(f"Number of trucks unloaded in an {SIM_TIME / 60:.0f} hour day = {unloaded_trucks_count[0]}")

# Execute the simulation and print the result.
solve_logistics_simulation()
<<<50>>>