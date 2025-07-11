import simpy
import random

def run_simulation():
    """
    This function sets up and runs the DC simulation.
    """
    # --- 1. Simulation Parameters ---
    SIM_TIME = 8 * 60      # Simulation time in minutes (8 hours)
    NUM_DOCKS = 3          # Number of docks
    MEAN_ARRIVAL = 7       # Mean time between truck arrivals (minutes)
    MIN_SERVICE = 20       # Minimum service time (minutes)
    MAX_SERVICE = 35       # Maximum service time (minutes)
    TRAVEL_TIME = 1        # Travel time from yard to dock (minutes)
    RANDOM_SEED = 42       # Seed for reproducibility

    # --- Global counter for unloaded trucks ---
    unloaded_trucks_count = 0

    def truck(env, name, docks):
        """
        Models the process for a single truck: arrival, queuing, travel, service, and departure.
        """
        nonlocal unloaded_trucks_count
        
        # A truck requests a dock. It will wait here if all docks are busy.
        with docks.request() as request:
            yield request

            # Once a dock is granted, the truck travels for 1 minute.
            yield env.timeout(TRAVEL_TIME)

            # The truck is then serviced (unloaded).
            service_time = random.uniform(MIN_SERVICE, MAX_SERVICE)
            yield env.timeout(service_time)
            
            # Check if the truck finished unloading within the 8-hour workday.
            # env.now gives the current simulation time.
            if env.now <= SIM_TIME:
                unloaded_trucks_count += 1

    def source(env, docks):
        """
        Generates new trucks arriving at the DC.
        """
        truck_id = 0
        while True:
            # Wait for the next truck to arrive.
            # The time between arrivals is exponentially distributed.
            interarrival_time = random.expovariate(1.0 / MEAN_ARRIVAL)
            yield env.timeout(interarrival_time)
            
            truck_id += 1
            # Start the process for the new truck.
            env.process(truck(env, f'Truck {truck_id}', docks))

    # --- 5. Run and Report ---
    random.seed(RANDOM_SEED)
    env = simpy.Environment()
    docks = simpy.Resource(env, capacity=NUM_DOCKS)
    env.process(source(env, docks))
    
    # Run the simulation for the 8-hour day.
    env.run(until=SIM_TIME)

    # Print the final result in the requested format.
    print(f"Total number of trucks unloaded in an 8-hour day = {unloaded_trucks_count}")
    return unloaded_trucks_count

# Execute the simulation and capture the result for the final answer.
final_count = run_simulation()
print(f"<<<{final_count}>>>")