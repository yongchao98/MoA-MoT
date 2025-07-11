import sys
import io

# Redirect stdout to capture print output for the final answer
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve():
    """
    Analyzes different insect behavioral patterns to determine which has the
    greatest positive effect on milkweed plant fitness.
    """
    print("Analyzing which behavioral pattern most benefits plant fitness.")
    print("----------------------------------------------------------------")
    print("Premise: Plant fitness is directly proportional to pollination success.")
    print("For milkweed, pollination occurs when an insect is feeding on nectar (behavior 5-6).")
    print("Therefore, the pattern that maximizes total feeding time has the greatest positive effect.\n")

    print("Interpreting the notation:")
    print(" - '6-5': Duration of a feeding bout.")
    print(" - '4-3': Duration of an interaction bout.")
    print(" - '2-1': Duration of an investigation bout.")
    print(" - An interaction (4-3) consists of feeding time (6-5) and non-feeding contact time.")
    print(" - The choice '6-5 >> 4-3' is interpreted as 'feeding time >> non-feeding interaction time'.\n")

    # --- Scenario B: Feeding Time (6-5) >> Non-Feeding Interaction Time ---
    # This insect is an efficient feeder, spending most of its contact time feeding.
    t_feed_B = 100
    t_non_feed_interact_B = 10
    t_interact_B = t_feed_B + t_non_feed_interact_B
    fitness_B = t_feed_B
    print("--- Pattern B: Feeding duration (6-5) >> Non-feeding interaction duration ---")
    print("This describes an efficient pollinator.")
    print(f"Equation: Fitness = Feeding Duration")
    print(f"Example values: Feeding duration = {t_feed_B}, Non-feeding interaction = {t_non_feed_interact_B}")
    print(f"Resulting Fitness Score = {fitness_B}\n")

    # --- Scenario A: Interaction Time (4-3) >> Feeding Time (6-5) ---
    # This insect is an inefficient feeder, spending most of its contact time not feeding.
    t_feed_A = 10
    t_non_feed_interact_A = 100
    t_interact_A = t_feed_A + t_non_feed_interact_A
    fitness_A = t_feed_A
    print("--- Pattern A: Interaction duration (4-3) >> Feeding duration (6-5) ---")
    print("This describes an inefficient pollinator.")
    print(f"Equation: Fitness = Feeding Duration")
    print(f"Example values: Feeding duration = {t_feed_A}, Non-feeding interaction = {t_non_feed_interact_A}")
    print(f"Resulting Fitness Score = {fitness_A}\n")

    # --- Scenario C: Interaction Time (4-3) >> Investigation Time (2-1) ---
    # This insect gets to business quickly, but the quality of the interaction is unknown.
    t_investigate_C = 10
    t_interact_C = 110  # Much greater than t_investigate_C
    # Assume a neutral 50/50 split for interaction time.
    t_feed_C = t_interact_C / 2
    fitness_C = t_feed_C
    print("--- Pattern C: Interaction duration (4-3) >> Investigation duration (2-1) ---")
    print("This insect spends lots of time on the plant, but may or may not be feeding.")
    print(f"Equation: Fitness = Feeding Duration (assuming 50% of interaction is feeding)")
    print(f"Example values: Total interaction = {t_interact_C}, Investigation = {t_investigate_C}")
    print(f"Resulting Fitness Score = {t_interact_C} / 2 = {fitness_C}\n")
    
    # --- Scenario E: n(1)/hour >> n(3)/hour ---
    # This insect is a "window shopper" that rarely makes contact.
    n1_E = 100  # investigations per hour
    n3_E = 10   # interactions per hour
    # Use the average feeding time from our efficient insect in B for comparison
    avg_feed_per_interaction_E = t_feed_B 
    total_feed_time_E = n3_E * avg_feed_per_interaction_E
    fitness_E = total_feed_time_E
    print("--- Pattern E: n(investigations) >> n(interactions) ---")
    print("This insect rarely lands on the plant.")
    print(f"Equation: Fitness = Interactions per hour * Avg. Feeding Duration per Interaction")
    print(f"Example values: Interactions per hour = {n3_E}, Avg. Feeding Duration = {avg_feed_per_interaction_E}")
    print(f"Resulting Fitness Score = {n3_E} * {avg_feed_per_interaction_E} = {fitness_E}\n")
    # Note: This score is for an hour, while others are per-event. But even if the insect from B only visits 10 times, its fitness score would be 10 * 100 = 1000, still higher than E.
    # Let's normalize E to a single event for better comparison.
    fitness_E_per_event = avg_feed_per_interaction_E * (n3_E / n1_E)
    print(f"For a fairer comparison, the expected fitness per investigation event is low: {fitness_E_per_event:.2f}\n")


    print("----------------------------------------------------------------")
    print("Conclusion:")
    print(f"Comparing the fitness scores per event: A={fitness_A}, B={fitness_B}, C={fitness_C}, E={fitness_E_per_event:.2f}")
    print("Pattern B yields the highest fitness score because it describes an insect that maximizes the time spent on the most crucial behavior for pollination: feeding.")
    print("While other patterns like C are positive, they do not guarantee efficient pollination in the same way B does.")

solve()
# Get the captured output and print it to the actual stdout
output = captured_output.getvalue()
sys.stdout = old_stdout
print(output)
print("<<<B>>>")