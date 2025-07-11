import math
from itertools import chain, combinations

def get_required_u_for_profile(t, m, agent_groups):
    """
    Calculates the minimum u required for a given preference profile.

    Args:
        t (int): The threshold for Condition 1.
        m (int): The number of items.
        agent_groups (dict): A dictionary where keys are preference orders (tuples)
                             and values are the number of agents with that preference.
    """
    items = list(range(1, m + 1))
    
    # Generate all possible subsets O of the items
    all_subsets = list(chain.from_iterable(combinations(items, r) for r in range(len(items) + 1)))

    min_required_u_for_profile = float('inf')

    print(f"Analyzing profile with t={t}, m={m}")
    print("Agent Groups:")
    for pref, count in agent_groups.items():
        print(f"  - {count} agents with preference {pref}")
    print("-" * 30)

    # Analyze each subset O
    for O in all_subsets:
        O = set(O)
        
        # --- Check Condition 1 ---
        cond1_satisfied = True
        if not O: # Condition 1 is vacuously true for the empty set
            pass
        else:
            for j in O:
                # Count agents whose favorite item in O is j
                count_j_in_O = 0
                for pref, num_agents in agent_groups.items():
                    # Find the first item from O in the preference list
                    favorite_in_O = next((item for item in pref if item in O), None)
                    if favorite_in_O == j:
                        count_j_in_O += num_agents
                
                if count_j_in_O <= t:
                    cond1_satisfied = False
                    break
        
        # If Condition 1 fails, this O can never be suitable.
        if not cond1_satisfied:
            print(f"O = {O if O else '{}'}: Fails Condition 1. Cannot be suitable.")
            continue

        # --- If Condition 1 holds, calculate u needed for Condition 2 ---
        max_prefer_k = 0
        items_not_in_O = [k for k in items if k not in O]

        if not items_not_in_O: # Condition 2 is vacuously true if O is the set of all items
            pass
        else:
            for k in items_not_in_O:
                # Count agents who prefer k over all items in O
                prefer_k_over_O = 0
                for pref, num_agents in agent_groups.items():
                    # Check if k appears before all items of O in the preference list
                    is_k_preferred = True
                    if not O: # k is preferred over an empty set of items
                         # This case corresponds to agents whose top choice is k
                        if pref[0] == k:
                            prefer_k_over_O += num_agents
                        continue

                    k_pos = pref.index(k)
                    for item_in_O in O:
                        if pref.index(item_in_O) < k_pos:
                            is_k_preferred = False
                            break
                    if is_k_preferred:
                        prefer_k_over_O += num_agents
                
                if prefer_k_over_O > max_prefer_k:
                    max_prefer_k = prefer_k_over_O
        
        required_u_for_O = max_prefer_k
        print(f"O = {O if O else '{}'}: Condition 1 OK. Requires u >= {required_u_for_O}.")
        
        if required_u_for_O < min_required_u_for_profile:
            min_required_u_for_profile = required_u_for_O

    print("-" * 30)
    print(f"For this profile, a suitable O can be found if u is at least {min_required_u_for_profile}.")
    print("This implies the answer to the problem must be at least this value.")
    print(f"The final answer is {min_required_u_for_profile}.")


if __name__ == '__main__':
    t_param = 20
    m_param = 4
    
    # The constructed "worst-case" preference profile
    worst_case_profile = {
        (1, 2, 3, 4): 20,
        (2, 3, 1, 4): 20,
        (3, 1, 2, 4): 20,
        (4, 1, 2, 3): 40,
    }
    
    get_required_u_for_profile(t_param, m_param, worst_case_profile)