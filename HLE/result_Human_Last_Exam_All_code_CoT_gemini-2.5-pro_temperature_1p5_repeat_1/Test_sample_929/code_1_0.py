import math

def calculate_fitness(description, total_time, a_rate, b_rate, a_duration, b_duration):
    """
    Calculates fitness based on relative rates of two behaviors.
    Fitness is measured by the number of feeding initiations, n(5).
    """
    print(f"Analyzing Pattern: {description}")
    
    # The ratio of behavior A to behavior B
    ratio = a_rate / b_rate
    
    # Total number of bouts that can fit in the total time
    # This is Total Time / (weighted average duration of a bout)
    total_bouts = total_time / ( (ratio * a_duration + 1 * b_duration) / (ratio + 1) )
    
    # Number of bouts of each type
    num_a_bouts = total_bouts * (ratio / (ratio + 1))
    num_b_bouts = total_bouts * (1 / (ratio + 1))

    # Pattern D is about the rate of feeding (n(5)) vs interaction (n(3))
    # Fitness is directly represented by n(5)
    fitness_score = num_a_bouts
    
    print(f"Description: The rate of Behavior A is {ratio:.1f} times the rate of Behavior B.")
    print(f"In {total_time} seconds, this results in:")
    print(f" - Number of Behavior A Bouts (e.g., n(5) feeding): {math.floor(num_a_bouts)}")
    print(f" - Number of Behavior B Bouts (e.g., n(3) interaction): {math.floor(num_b_bouts)}")
    print(f"Fitness Score (proportional to n(5)): {math.floor(fitness_score)}\n")
    return fitness_score

def calculate_fitness_duration(description, total_time, a_duration_share, b_duration_share, a_bout_duration, b_bout_duration):
    """
    Calculates fitness based on relative durations of two behaviors.
    """
    print(f"Analyzing Pattern: {description}")
    
    # Total time spent on each behavior
    total_a_time = total_time * (a_duration_share / (a_duration_share + b_duration_share))
    total_b_time = total_time * (b_duration_share / (a_duration_share + b_duration_share))

    # Number of bouts of behavior A
    num_a_bouts = total_a_time / a_bout_duration
    
    # In this model, fitness is n(5), which is the number of feeding bouts.
    # Behavior A represents feeding.
    fitness_score = num_a_bouts

    print(f"Description: The duration of Behavior A is {a_duration_share / b_duration_share:.1f} times the duration of Behavior B.")
    print(f"In {total_time} seconds, this results in:")
    print(f" - Time on Behavior A (e.g., feeding): {total_a_time:.0f}s")
    print(f" - Time on Behavior B (e.g., interaction): {total_b_time:.0f}s")
    print(f"Fitness Score (n(5) from feeding bouts): {math.floor(fitness_score)}\n")
    return fitness_score


# --- Simulation Parameters ---
VISIT_TIME_SECONDS = 100  # Total time the insect is in contact with the plant
FEEDING_BOUT_DURATION = 8  # 8 seconds to feed from one flower
INTERACTION_BOUT_DURATION = 2 # 2 seconds to walk between flowers

print("--- Pollination Fitness Analysis ---")
print(f"This simulation assumes a total contact time of {VISIT_TIME_SECONDS}s.")
print("Plant fitness is determined by the number of flowers an insect visits for feeding.")
print("Therefore, the highest Fitness Score (number of feeding bouts, n(5)) is the best outcome for the plant.\n")

# --- Evaluate Each Plausible Choice ---

# Choice A: 4-3 >> 6-5 (duration of interaction >> duration of feeding)
# Let interaction duration be 9x feeding duration
# Behavior A is feeding, Behavior B is interaction. We are swapping the duration shares.
score_A = calculate_fitness_duration(
    "A. dur(interaction) >> dur(feeding)", VISIT_TIME_SECONDS, 
    a_duration_share=1, b_duration_share=9,
    a_bout_duration=FEEDING_BOUT_DURATION, 
    b_bout_duration=INTERACTION_BOUT_DURATION
)

# Choice B: 6-5 >> 4-3 (duration of feeding >> duration of interaction)
# Let feeding duration be 9x interaction duration
score_B = calculate_fitness_duration(
    "B. dur(feeding) >> dur(interaction)", VISIT_TIME_SECONDS, 
    a_duration_share=9, b_duration_share=1,
    a_bout_duration=FEEDING_BOUT_DURATION, 
    b_bout_duration=INTERACTION_BOUT_DURATION
)

# Choice D: n(5)/hour >> n(3)/hour (rate of feeding >> rate of interaction)
# Let the rate of feeding bouts be 9x the rate of interaction bouts
score_D = calculate_fitness(
    "D. n(feeding)/hr >> n(interaction)/hr", VISIT_TIME_SECONDS,
    a_rate=9, b_rate=1, 
    a_duration=FEEDING_BOUT_DURATION, 
    b_duration=INTERACTION_BOUT_DURATION
)

print("--- Conclusion ---")
print("Choices C, E, and F are not modeled directly:")
print("- C (dur(interaction) >> dur(investigation)) is good, but less specific about pollination than B or D.")
print("- E (n(investigation) >> n(interaction)) is bad, as it means few contacts are made.")
print("- F (n(interaction) >> n(investigation)) is logically impossible.")

print("\nComparing the fitness scores from plausible scenarios:")
print(f"Score for A: {math.floor(score_A)}")
print(f"Score for B: {math.floor(score_B)}")
print(f"Score for D: {math.floor(score_D)}")

if score_D > score_B and score_D > score_A:
    print("\nPattern D yields the highest fitness score. A high frequency of feeding initiations (n(5)) means the insect visits many different flowers, maximizing the chances of cross-pollination.")
    final_answer = 'D'
else:
    # This case shouldn't be reached with the current model, but as a fallback:
    print("\nBased on the simulation, the best pattern needs to be re-evaluated.")
    final_answer = 'Error'
    
# The final result as requested by the user format
# Note: The code is just a model to prove the reasoning. The reasoning itself is the core of the solution.
# The reasoning concludes that maximizing the number of flower visits (n(5)) is best.
# Pattern D describes a high frequency of feeding initiations (n(5)) compared to other contact (n(3)).
# This behavior directly translates to visiting many flowers.
# Thus, Pattern D has the greatest positive effect on plant fitness.

print("\nThe equation representing the winning pattern is: n(5)/hour >> n(3)/hour")

<<<D>>>