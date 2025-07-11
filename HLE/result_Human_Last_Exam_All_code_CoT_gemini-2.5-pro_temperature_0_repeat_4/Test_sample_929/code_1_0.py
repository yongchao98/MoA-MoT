import collections

def calculate_metrics(events, total_duration):
    """
    Calculates total investigation, interaction, and feeding times from a list of events.
    An event is a tuple of (time, code).
    Ethogram:
    1: investigation start, 2: investigation end
    3: interaction start, 4: interaction end
    5: feeding start, 6: feeding end
    """
    # Helper to find matching start/end events
    def get_total_duration(start_code, end_code, events):
        total_time = 0
        starts = sorted([t for t, c in events if c == start_code])
        ends = sorted([t for t, c in events if c == end_code])
        
        # Simple 1-to-1 matching for this model
        for start_time, end_time in zip(starts, ends):
            if end_time > start_time:
                total_time += (end_time - start_time)
        return total_time

    # Handle implied events
    # 3 after 1 implies 2: An interaction starting implies an investigation ending.
    # We model this by making the start of interaction the end of investigation.
    # 4 after 5 implies 6: An interaction ending implies feeding ends.
    # We model this by making the end of interaction the end of feeding.
    
    # For simplicity in this model, we'll use explicit start/end pairs in our data.

    t_investigation = get_total_duration(1, 2, events)
    t_interaction = get_total_duration(3, 4, events)
    t_feeding = get_total_duration(5, 6, events)
    
    # Frequencies
    n_investigation = len([c for t, c in events if c == 1])
    n_interaction = len([c for t, c in events if c == 3])
    n_feeding = len([c for t, c in events if c == 5])

    # Fitness score is the total feeding time
    fitness_score = t_feeding
    
    return {
        "T_investigation": t_investigation,
        "T_interaction": t_interaction,
        "T_feeding": t_feeding,
        "N_investigation": n_investigation,
        "N_interaction": n_interaction,
        "N_feeding": n_feeding,
        "Fitness_Score": fitness_score
    }

def run_simulation():
    """
    Simulates different behavioral patterns and calculates their fitness effect.
    """
    scenarios = collections.OrderedDict()
    
    # Scenario A: 4-3 >> 6-5 (Long interaction, short feeding)
    # Insect loiters on the plant but doesn't feed much. Low efficiency.
    scenarios['A. Interaction >> Feeding'] = [
        (0, 1), (10, 3),              # Investigation: 10s
        (20, 5), (30, 6),             # Feeding: 10s
        (200, 4), (210, 2)            # Interaction: 190s
    ]

    # Scenario B: 6-5 >> 4-3 (Feeding duration is the majority of interaction duration)
    # Highly efficient pollinator. Gets on the plant and feeds the whole time.
    scenarios['B. Feeding ≈ Interaction'] = [
        (0, 1), (10, 3),              # Investigation: 10s
        (12, 5), (198, 6),            # Feeding: 186s
        (200, 4), (210, 2)            # Interaction: 190s
    ]

    # Scenario C: 4-3 >> 2-1 (Long interaction, short investigation)
    # Insect commits to a long visit once it arrives.
    scenarios['C. Interaction >> Investigation'] = [
        (0, 1), (5, 3),               # Investigation: 5s
        (20, 5), (120, 6),            # Feeding: 100s
        (200, 4), (205, 2)            # Interaction: 195s
    ]

    # Scenario E: n(1)/hr >> n(3)/hr (Many investigations, few interactions)
    # "Window shopper" insect. Hesitant to land.
    events_e = []
    for i in range(9): # 9 investigations with no contact
        events_e.extend([(i*20, 1), (i*20 + 10, 2)])
    # 1 investigation that leads to contact
    events_e.extend([(180, 1), (190, 3), (195, 5), (205, 6), (210, 4), (220, 2)])
    scenarios['E. N(investigate) >> N(interact)'] = events_e

    # Scenario F: n(3)/hr >> n(1)/hr (Many interactions per investigation)
    # Insect makes many short contacts during one visit to the plant area.
    events_f = [(0, 1)]
    for i in range(5): # 5 separate interactions
        base_time = 10 + i * 40
        events_f.extend([
            (base_time, 3), (base_time + 5, 5), 
            (base_time + 15, 6), (base_time + 20, 4)
        ])
    events_f.append((220, 2))
    scenarios['F. N(interact) >> N(investigate)'] = events_f

    print("--- Behavioral Pattern Fitness Analysis ---")
    print("Fitness Score is defined as total feeding time.\n")
    
    results = {}
    for name, events in scenarios.items():
        metrics = calculate_metrics(events, 3600)
        results[name] = metrics
        
        print(f"Scenario: {name}")
        print(f"  Total Interaction Time: {metrics['T_interaction']:>4}s")
        print(f"  Total Feeding Time:     {metrics['T_feeding']:>4}s")
        
        # Calculate efficiency ratio to illustrate the point
        efficiency = 0
        if metrics['T_interaction'] > 0:
            efficiency = (metrics['T_feeding'] / metrics['T_interaction']) * 100
        print(f"  Feeding Efficiency:     {efficiency:5.1f}%")
        print(f"  FITNESS SCORE:          {metrics['Fitness_Score']:>4}\n")

    # Find the best scenario
    best_scenario = max(results, key=lambda k: results[k]['Fitness_Score'])
    print("--- Conclusion ---")
    print(f"The pattern with the greatest positive effect on plant fitness is '{best_scenario}',")
    print("as it results in the highest total feeding time.")

run_simulation()
<<<B>>>