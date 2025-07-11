import pandas as pd

def calculate_pollination_fitness(scenario_name, n_investigate, p_interact, p_feed, d_feed, d_interact_total):
    """
    Calculates and prints a pollination fitness score based on behavioral parameters.
    
    Formula: Pollination Fitness Score = (Number of Feeding Bouts) * (Average Feeding Duration)
    
    Args:
        scenario_name (str): The name of the scenario (e.g., 'Baseline', 'A').
        n_investigate (int): Number of investigations per hour.
        p_interact (float): Probability of an investigation leading to an interaction (n(3)/n(1)).
        p_feed (float): Probability of an interaction leading to feeding (n(5)/n(3)).
        d_feed (float): Average duration of a feeding bout in seconds (6-5).
        d_interact_total (float): Average total duration of an interaction in seconds (4-3).
    """
    n_interact = n_investigate * p_interact
    n_feed = n_interact * p_feed
    
    pollination_score = n_feed * d_feed
    
    print(f"Scenario {scenario_name}:")
    print(f"  - Behavioral Pattern: {pattern_descriptions[scenario_name]}")
    print(f"  - Number of feeding bouts = {n_investigate} investigations * {p_interact:.2f} P(interact) * {p_feed:.2f} P(feed) = {n_feed:.1f} bouts")
    print(f"  - Average feeding duration = {d_feed:.1f} seconds")
    print(f"  - Pollination Fitness = {n_feed:.1f} feeding bouts * {d_feed:.1f} seconds/bout = {pollination_score:.1f}")
    print("-" * 20)
    
    return pollination_score

# --- Scenario Setup ---
# A researcher observes 100 investigations per hour as a baseline.
n_investigate_base = 100

pattern_descriptions = {
    'Baseline': 'A neutral or average scenario for comparison.',
    'A': '4-3 >> 6-5 (Interaction duration is much longer than feeding duration)',
    'B': '6-5 >> 4-3 (Feeding duration is much longer than other interaction time)',
    'C': '4-3 >> 2-1 (Interaction duration is much longer than investigation duration)',
    'D': 'n(5)/hr >> n(3)/hr (High probability of feeding given an interaction)',
    'E': 'n(1)/hr >> n(3)/hr (Low probability of interaction given an investigation)',
    'F': 'n(3)/hr >> n(1)/hr (High probability of interaction given an investigation)'
}

# --- Parameter Definitions for each scenario ---

# Baseline: Averages for all parameters
scenarios = {
    'Baseline': {'p_interact': 0.5, 'p_feed': 0.5, 'd_feed': 10.0, 'd_interact_total': 20.0},
    # A: Interaction is long, but feeding is very short.
    'A': {'p_interact': 0.5, 'p_feed': 0.5, 'd_feed': 2.0, 'd_interact_total': 20.0},
    # B: Interaction time is almost entirely feeding time. This is the highest quality visit.
    'B': {'p_interact': 0.5, 'p_feed': 0.5, 'd_feed': 18.0, 'd_interact_total': 20.0},
    # C: Interaction is long. We'll keep feeding as half of it, similar to baseline.
    'C': {'p_interact': 0.5, 'p_feed': 0.5, 'd_feed': 15.0, 'd_interact_total': 30.0},
    # D: High conversion from interaction to feeding.
    'D': {'p_interact': 0.5, 'p_feed': 0.9, 'd_feed': 10.0, 'd_interact_total': 20.0},
    # E: Low conversion from investigation to interaction (most visitors don't land).
    'E': {'p_interact': 0.1, 'p_feed': 0.5, 'd_feed': 10.0, 'd_interact_total': 20.0},
    # F: High conversion from investigation to interaction (most visitors land).
    'F': {'p_interact': 0.9, 'p_feed': 0.5, 'd_feed': 10.0, 'd_interact_total': 20.0}
}

# --- Calculation and Output ---
print("Calculating Pollination Fitness Score for Different Behavioral Scenarios.")
print("The greatest positive effect on plant fitness corresponds to the highest score.\n")
results = {}
for name, params in scenarios.items():
    score = calculate_pollination_fitness(
        scenario_name=name,
        n_investigate=n_investigate_base,
        p_interact=params['p_interact'],
        p_feed=params['p_feed'],
        d_feed=params['d_feed'],
        d_interact_total=params['d_interact_total']
    )
    results[name] = score

# --- Final Conclusion ---
best_scenario = max(results, key=results.get)
print(f"\nConclusion: Scenario '{best_scenario}' has the highest pollination fitness score ({results[best_scenario]:.1f}).")
print(f"This corresponds to the behavioral pattern: '{pattern_descriptions[best_scenario]}'.")
print("This pattern describes a highly efficient pollinator that spends the majority of its contact time actively feeding, which is the most beneficial action for the plant's reproduction.")
