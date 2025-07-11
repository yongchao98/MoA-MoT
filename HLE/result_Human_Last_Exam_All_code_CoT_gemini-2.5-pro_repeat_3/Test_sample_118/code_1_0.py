import itertools

def solve_coin_problem():
    """
    Analyzes a strategy for the coin problem to find the guaranteed number of real coins.

    The strategy is as follows:
    1. Partition 1000 coins into A(332), B(332), C(332), D(4).
    2. Weighing 1: A vs B.
    3. Weighing 2: B vs C.

    The function iterates through all 9 possible outcomes of the two weighings,
    determines the possible distributions of the 4 fake coins for each outcome,
    and calculates how many coins can be guaranteed as real. The final answer
    is the minimum number of guaranteed coins over all possible outcomes.
    """
    group_sizes = {'A': 332, 'B': 332, 'C': 332, 'D': 4}
    groups = list(group_sizes.keys())
    total_fakes = 4

    # The outcomes are represented by strings: 'gt' (>), 'lt' (<), 'eq' (=)
    outcomes_w1 = ['gt', 'lt', 'eq']  # A vs B
    outcomes_w2 = ['gt', 'lt', 'eq']  # B vs C
    
    min_guaranteed_coins = float('inf')

    print("Analyzing the strategy: Partition into A(332), B(332), C(332), D(4).")
    print("Weighing 1: A vs B. Weighing 2: B vs C.\n")

    # Iterate through all 9 possible outcomes
    for o1 in outcomes_w1:
        for o2 in outcomes_w2:
            
            possible_distributions = []
            
            # Generate all possible distributions of 4 fakes into 4 groups
            # This is equivalent to finding non-negative integer solutions to fA+fB+fC+fD=4
            for f_A in range(total_fakes + 1):
                for f_B in range(total_fakes - f_A + 1):
                    for f_C in range(total_fakes - f_A - f_B + 1):
                        f_D = total_fakes - f_A - f_B - f_C
                        
                        dist = {'A': f_A, 'B': f_B, 'C': f_C, 'D': f_D}
                        
                        # Check if distribution is possible given group sizes
                        if any(dist[g] > group_sizes[g] for g in groups):
                            continue

                        # Check if distribution matches the weighing outcomes
                        w1_match = False
                        if o1 == 'gt' and dist['A'] < dist['B']: w1_match = True
                        if o1 == 'lt' and dist['A'] > dist['B']: w1_match = True
                        if o1 == 'eq' and dist['A'] == dist['B']: w1_match = True

                        w2_match = False
                        if o2 == 'gt' and dist['B'] < dist['C']: w2_match = True
                        if o2 == 'lt' and dist['B'] > dist['C']: w2_match = True
                        if o2 == 'eq' and dist['B'] == dist['C']: w2_match = True

                        if w1_match and w2_match:
                            possible_distributions.append(dist)

            # Determine the number of guaranteed real coins for this outcome
            guaranteed_coins_this_outcome = 0
            if possible_distributions:
                # A group is guaranteed real if it contains 0 fakes in ALL possible distributions
                for group_name in groups:
                    is_guaranteed_real = True
                    for dist in possible_distributions:
                        if dist[group_name] != 0:
                            is_guaranteed_real = False
                            break
                    if is_guaranteed_real:
                        guaranteed_coins_this_outcome += group_sizes[group_name]

            # Update the minimum guarantee
            min_guaranteed_coins = min(min_guaranteed_coins, guaranteed_coins_this_outcome)

            # Print the analysis for the current outcome
            outcome_str = f"Outcome: A {'>' if o1=='gt' else '<' if o1=='lt' else '='} B, B {'>' if o2=='gt' else '<' if o2=='lt' else '='} C"
            print(f"{outcome_str}")
            if not possible_distributions:
                # This can happen if an outcome is impossible (e.g., A>B and B>C and C>A)
                # An impossible outcome provides no information but doesn't limit our guarantee.
                # We can consider the guarantee infinite for this branch, so it doesn't affect the minimum.
                print("  This outcome is impossible.")
                print(f"  Guaranteed real coins for this branch: N/A")
            else:
                # Find which groups are identified as real
                identified_groups = []
                if guaranteed_coins_this_outcome > 0:
                    for group_name in groups:
                        if all(dist[group_name] == 0 for dist in possible_distributions):
                            identified_groups.append(f"{group_name}({group_sizes[group_name]})")
                
                if identified_groups:
                    print(f"  Unique state found or group(s) {', '.join(identified_groups)} are always fake-free.")
                else:
                    print(f"  Ambiguity in fake coin locations prevents a guarantee.")
                
                print(f"  Guaranteed real coins for this branch: {guaranteed_coins_this_outcome}")
            print("-" * 20)


    print(f"\nThe maximum number of real coins we can GUARANTEE to identify is the minimum of the guarantees from all possible outcomes.")
    print(f"Final Answer: The minimum guarantee across all 9 outcomes is {min_guaranteed_coins}.")
    
solve_coin_problem()