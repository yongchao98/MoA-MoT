import pandas as pd

def simulate_pollination(pattern):
    """
    Simulates insect behavior based on a pattern and calculates a pollination score.
    """
    # Plant setup: an umbel with 10 flowers
    num_flowers = 10
    flowers = {i: {'has_pollen': True, 'is_pollinated': False} for i in range(num_flowers)}
    
    # Insect state
    insect_pollen_carrier = 0
    pollination_score = 0
    
    # Ethogram counters
    n_investigate = 0
    n_interact = 0
    n_feed = 0
    t_investigate = 0
    t_interact = 0
    t_feed = 0

    # The core logic of pollination: moving between flowers
    def visit_flower(flower_id):
        nonlocal insect_pollen_carrier, pollination_score
        # Deposit pollen if carrying any and flower is not yet pollinated
        if insect_pollen_carrier > 0 and not flowers[flower_id]['is_pollinated']:
            flowers[flower_id]['is_pollinated'] = True
            pollination_score += 1
            insect_pollen_carrier -= 1 # Pollen is deposited
        # Pick up new pollen
        if flowers[flower_id]['has_pollen']:
            insect_pollen_carrier += 1
            flowers[flower_id]['has_pollen'] = False

    # --- Simulation based on pattern ---
    
    if pattern == 'A': # t_interact >> t_feed
        n_investigate = 1
        t_investigate = 10
        n_interact = 1
        t_interact = 100
        n_feed = 1
        t_feed = 5
        # Long interaction on one flower with short feeding
        visit_flower(0)

    elif pattern == 'C': # t_interact >> t_investigate
        n_investigate = 1
        t_investigate = 5
        n_interact = 1
        t_interact = 100
        n_feed = 1
        t_feed = 50 # Assumes long interaction includes significant feeding
        # Long interaction on one flower
        visit_flower(0)

    elif pattern == 'E': # n_investigate >> n_interact
        n_investigate = 10
        t_investigate = 10 * 10
        # Only one investigation leads to a brief interaction
        n_interact = 1
        t_interact = 10
        n_feed = 1
        t_feed = 5
        visit_flower(0)

    elif pattern == 'F': # n_interact >> n_investigate
        n_investigate = 1
        t_investigate = 100 # One long visit period
        # This one visit has many interactions
        n_interact = 10 
        for i in range(n_interact):
            t_interact_bout = 10
            t_feed_bout = 8
            
            t_interact += t_interact_bout
            t_feed += t_feed_bout
            n_feed += 1
            
            # Each interaction is on a new flower
            visit_flower(i % num_flowers)
            
    # For impossible scenarios B and D, we show zero score
    elif pattern in ['B', 'D']:
        pass


    # Format results nicely for printing
    results = {
        "Pattern": pattern,
        "n(investigate)": n_investigate,
        "t(investigate)": t_investigate,
        "n(interact)": n_interact,
        "t(interact)": t_interact,
        "n(feed)": n_feed,
        "t(feed)": t_feed,
        "Pollination Score": pollination_score
    }
    return results

def run_analysis():
    """
    Runs the simulation for all patterns and prints the results.
    """
    patterns = ['A', 'B', 'C', 'D', 'E', 'F']
    all_results = []

    print("Simulating behavioral patterns to assess impact on plant fitness (pollination)...")
    print("-" * 70)

    for p in patterns:
        res = simulate_pollination(p)
        
        # Build readable description of the scenario
        desc = ""
        if p == 'A': desc = "4-3 >> 6-5 (Long interaction, short feeding)"
        elif p == 'B': desc = "6-5 >> 4-3 (Impossible: t_feed > t_interact)"
        elif p == 'C': desc = "4-3 >> 2-1 (Long interaction, short investigation)"
        elif p == 'D': desc = "n(5)/hr >> n(3)/hr (Impossible: n_feed > n_interact)"
        elif p == 'E': desc = "n(1)/hr >> n(3)/hr (Many investigations, few interactions)"
        elif p == 'F': desc = "n(3)/hr >> n(1)/hr (Many interactions per investigation)"

        print(f"Scenario {p}: {desc}")
        
        # Printing the numbers in the "equation"
        if p == 'A':
            print(f"t(interaction) = {res['t(interact)']}, t(feeding) = {res['t(feed)']}. Check: {res['t(interact)']} >> {res['t(feed)']}")
        elif p == 'C':
            print(f"t(interaction) = {res['t(interact)']}, t(investigation) = {res['t(investigate)']}. Check: {res['t(interact)']} >> {res['t(investigate)']}")
        elif p == 'E':
            print(f"n(investigation) = {res['n(investigate)']}, n(interaction) = {res['n(interact)']}. Check: {res['n(investigate)']} >> {res['n(interact)']}")
        elif p == 'F':
            print(f"n(interaction) = {res['n(interact)']}, n(investigation) = {res['n(investigate)']}. Check: {res['n(interact)']} >> {res['n(investigate)']}")
        
        print(f"Resulting Pollination Score: {res['Pollination Score']}\n")
        all_results.append(res)
    
    print("-" * 70)
    df = pd.DataFrame(all_results)
    print("Summary Table:")
    print(df.to_string(index=False))
    
    # Final conclusion based on the simulation
    max_score = df['Pollination Score'].max()
    best_pattern = df[df['Pollination Score'] == max_score]['Pattern'].iloc[0]
    print(f"\nConclusion: Pattern '{best_pattern}' results in the highest pollination score, indicating the greatest positive effect on plant fitness.")

if __name__ == '__main__':
    run_analysis()