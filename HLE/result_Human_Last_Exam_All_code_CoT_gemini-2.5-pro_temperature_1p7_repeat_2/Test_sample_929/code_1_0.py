def calculate_metrics(events):
    """
    Calculates total feeding time and other relevant metrics from a list of events.
    An event is a tuple of (time, event_code).
    """
    total_feeding_time = 0
    feeding_start_time = None
    
    # Store all times for each event code
    times = {code: [] for code in range(1, 7)}
    for time, code in sorted(events):
        times[code].append(time)
        if code == 5: # feeding start
            feeding_start_time = time
        elif code == 6: # feeding end
            if feeding_start_time is not None:
                duration = time - feeding_start_time
                total_feeding_time += duration
                feeding_start_time = None # Reset for next feeding bout

    # Calculate specific durations for comparison (assuming one bout per scenario for simplicity)
    investigation_duration = (times[2][0] - times[1][0]) if (times[1] and times[2]) else 0
    interaction_duration = (times[4][0] - times[3][0]) if (times[3] and times[4]) else 0
    # Use the first full feeding bout for comparison
    feeding_duration = (times[6][0] - times[5][0]) if (times[5] and times[6]) else 0
    
    return {
        'fitness_score': total_feeding_time,
        'investigation_duration': investigation_duration,
        'interaction_duration': interaction_duration,
        'feeding_duration': feeding_duration,
        'n1': len(times[1]),
        'n3': len(times[3]),
        'times': times
    }

def main():
    """
    Main function to evaluate all scenarios and determine the best one for plant fitness.
    """
    print("Analyzing behavioral patterns to find the one with the greatest positive effect on plant fitness.")
    print("Assumption: Plant fitness is directly proportional to pollination, which is best represented by the total time the insect spends feeding.\n")

    scenarios = {
        'A': {
            'description': "Interaction duration (4-3) >> Feeding duration (6-5)",
            'events': [
                (0, 1), (2, 3), # Investigation -> Interaction
                (10, 5), (15, 6), # Short feeding bout
                (60, 4) # End long interaction
            ]
        },
        'B': {
            'description': "Feeding duration (6-5) >> Interaction duration (4-3)",
            'events': [
                (0, 1), (2, 3), # Investigation -> Interaction
                (3, 5), # Start feeding almost immediately
                (58, 6), # End long feeding bout
                (60, 4) # End interaction shortly after
            ]
        },
        'C': {
            'description': "Interaction duration (4-3) >> Investigation duration (2-1)",
            'events': [
                (0, 1), (2, 2), # Short investigation
                (2, 3), # Immediately start interaction
                (10, 5), (25, 6), # Moderate feeding
                (60, 4) # End long interaction
            ]
        },
        'E': {
            'description': "n(1)/hour >> n(3)/hour",
            'events': [
                (0, 1), (5, 2), # Investigation 1, no interaction
                (10, 1), (15, 2), # Investigation 2, no interaction
                (20, 1), (25, 2), # Investigation 3, no interaction
                (30, 1), (32, 3), # Investigation 4 leads to interaction
                (35, 5), (40, 6), # Short feeding
                (45, 4) # End interaction
            ]
        }
    }

    results = {}

    for choice, scenario in scenarios.items():
        print(f"--- Scenario for Choice {choice}: {scenario['description']} ---")
        metrics = calculate_metrics(scenario['events'])
        results[choice] = metrics['fitness_score']
        
        if choice == 'A':
            t = metrics['times']
            print(f"Equation: ({t[4][0]} - {t[3][0]}) >> ({t[6][0]} - {t[5][0]})")
            print(f"Result: {metrics['interaction_duration']} >> {metrics['feeding_duration']}. The scenario holds.")
        
        elif choice == 'B':
            # This choice is interpreted as "most of interaction time is spent feeding"
            t = metrics['times']
            print("Interpretation: The time spent feeding makes up the vast majority of the time spent interacting.")
            print(f"Equation for feeding duration: {t[6][0]} - {t[5][0]} = {metrics['feeding_duration']}")
            print(f"Equation for interaction duration: {t[4][0]} - {t[3][0]} = {metrics['interaction_duration']}")
            print(f"Result: Feeding duration of {metrics['feeding_duration']} is the dominant part of the interaction duration of {metrics['interaction_duration']}. The scenario holds.")
        
        elif choice == 'C':
            t = metrics['times']
            print(f"Equation: ({t[4][0]} - {t[3][0]}) >> ({t[2][0]} - {t[1][0]})")
            print(f"Result: {metrics['interaction_duration']} >> {metrics['investigation_duration']}. The scenario holds.")

        elif choice == 'E':
            t = metrics['times']
            print(f"Equation: n(1) >> n(3)")
            print(f"Result: {metrics['n1']} investigations >> {metrics['n3']} interactions. The scenario holds.")
        
        print(f"Fitness Score (Total Feeding Time): {metrics['fitness_score']}\n")

    print("--- Impossible Scenarios ---")
    print("Choice D (n(5)/hour >> n(3)/hour) and F (n(3)/hour >> n(1)/hour) are logically impossible.")
    print("A feeding event (5) requires an interaction (3), and an interaction (3) requires an investigation (1).")
    print("Therefore, n(5) <= n(3) <= n(1) must always be true.\n")
    
    # Find the best choice
    best_choice = max(results, key=results.get)
    print("--- Conclusion ---")
    print("Comparing the fitness scores from the plausible scenarios:")
    for choice, score in results.items():
        print(f"  - Choice {choice}: {score}")
    
    print(f"\nScenario '{best_choice}' yields the highest fitness score. This behavioral pattern, where the insect spends the vast majority of its interaction time actively feeding, is the most beneficial for the plant's reproductive success.")

if __name__ == '__main__':
    main()
<<<B>>>