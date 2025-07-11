import sys
import io

# A simple class to structure the simulation and keep track of the state.
class PondSimulation:
    def __init__(self):
        self.time = 0
        self.goldfish_in_pond = 10
        self.goldfish_in_basket = 0
        self.shark_eaten_count = 0
        # We'll use a list to store event strings for the final log.
        self.event_log = []

    def log_event(self, timestamp, message):
        """Logs an event at a specific timestamp, updating the simulation time."""
        self.time = timestamp
        log_entry = f"Time: {self.time:2d} min | Pond GF: {self.goldfish_in_pond:2d} | Basket: {self.goldfish_in_basket} | Shark Eaten: {self.shark_eaten_count} | {message}"
        self.event_log.append(log_entry)

    def run_optimal_strategy(self):
        """Executes the step-by-step optimal strategy."""

        self.log_event(0, "Start. Fisherman uses Rod B (fastest catch).")

        # --- Phase 1: Initial Catch with Rod B (t=0 to t=15) ---
        # Goal: Catch fish quickly without filling the basket.
        
        # t=5: Fisherman catches his first fish.
        self.goldfish_in_pond -= 1
        self.goldfish_in_basket += 1
        self.log_event(5, "Fisherman catches 1st goldfish (Rod B).")

        # t=10: Two events happen simultaneously.
        self.goldfish_in_pond -= 1
        self.shark_eaten_count += 1
        self.log_event(10, "Shark eats 1st free-swimming goldfish.")
        
        self.goldfish_in_pond -= 1
        self.goldfish_in_basket += 1
        self.log_event(10, "Fisherman catches 2nd goldfish (Rod B).")

        # t=15: Fisherman catches his third fish. Basket now has 3.
        self.goldfish_in_pond -= 1
        self.goldfish_in_basket += 1
        self.log_event(15, "Fisherman catches 3rd goldfish (Rod B). Basket at 3.")
        phase1_duration = 15

        # --- Phase 2: Power-up the Shark (t=15 to t=24) ---
        # Goal: Use Rod A to feed the shark, empty the basket, and get the shark to "fast mode".
        self.log_event(15, "Switches to Rod A to feed the shark from the basket.")

        # Fisherman feeds the 3 fish from the basket. Each feeding takes 2 minutes.
        self.goldfish_in_basket -= 1
        self.shark_eaten_count += 1
        self.log_event(17, "Feeds shark with 1st fish from basket.")

        self.goldfish_in_basket -= 1
        self.shark_eaten_count += 1
        self.log_event(19, "Feeds shark with 2nd fish from basket.")

        # t=20: The shark's natural feeding timer is up.
        self.goldfish_in_pond -= 1
        self.shark_eaten_count += 1
        self.log_event(20, "Shark eats 2nd free-swimming goldfish.")

        # The third feeding finishes at t=21 (started at t=19).
        self.goldfish_in_basket -= 1
        self.shark_eaten_count += 1
        self.log_event(21, "Feeds shark with 3rd fish. Basket is empty.")

        # The fisherman has been using Rod A since t=15. A catch happens after 7 minutes.
        self.goldfish_in_pond -= 1
        self.goldfish_in_basket += 1
        self.log_event(22, "Catches 4th goldfish (Rod A). Rod B timer is reset.")

        # Feed this new fish to get the shark to its 6th meal.
        self.goldfish_in_basket -= 1
        self.shark_eaten_count += 1
        self.log_event(24, "Feeds shark with 4th fish. Shark has eaten 6 fish now!")
        phase2_duration = 24 - 15

        # --- Phase 3: Final Cleanup (t=24 to t=30) ---
        # Goal: Use the fast shark and fast Rod B to clear the remaining fish.
        self.log_event(24, "SHARK IS NOW FAST (eats every 2 min). Switches back to Rod B.")

        # The fast shark eats every 2 mins. Last ate at t=24.
        self.goldfish_in_pond -= 1
        self.shark_eaten_count += 1
        self.log_event(26, "Fast shark eats a goldfish.")

        self.goldfish_in_pond -= 1
        self.shark_eaten_count += 1
        self.log_event(28, "Fast shark eats another goldfish.")

        # Fisherman's Rod B catch, started at t=24, finishes after 5 minutes.
        self.goldfish_in_pond -= 1
        self.goldfish_in_basket += 1
        self.log_event(29, "Fisherman catches 5th goldfish (Rod B).")

        # The shark's final meal.
        self.goldfish_in_pond -= 1
        self.shark_eaten_count += 1
        self.log_event(30, "Fast shark eats the last goldfish in the pond.")
        phase3_duration = 30 - 24
        
        # --- Final Output ---
        print("The optimal strategy has been simulated. Here is the event log:\n")
        for entry in self.event_log:
            print(entry)
        
        final_time = self.time
        print("\n--------------------------------------------------------------")
        print("           FINAL CALCULATION")
        print("--------------------------------------------------------------")
        print("The earliest possible time when there are no more fish left is the sum of the three strategic phases.")
        print(f"Phase 1 Duration: {phase1_duration} minutes")
        print(f"Phase 2 Duration: {phase2_duration} minutes")
        print(f"Phase 3 Duration: {phase3_duration} minutes")
        print("\nFinal Equation:")
        print(f"{phase1_duration} + {phase2_duration} + {phase3_duration} = {final_time}")
        print(f"\nThe earliest time is {final_time} minutes.")


# Run the simulation and print the results.
simulation = PondSimulation()
simulation.run_optimal_strategy()