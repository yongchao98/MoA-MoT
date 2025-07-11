from itertools import combinations

def format_fluent(f_str):
    """Converts an internal string representation like 'f0(o1,o2)' to the required output format."""
    name_params = f_str.replace(')', '').split('(')
    name = name_params[0]
    params = name_params[1].split(',')
    
    fluent_map = {'f0': 'fluent_0', 'f1': 'fluent_1', 'f2': 'fluent_2'}
    obj_map = {'o1': 'object_1_type_0', 'o2': 'object_2_type_0'}
    
    formatted_name = fluent_map[name]
    formatted_params = [obj_map[p.strip()] for p in params]
    
    return f"{formatted_name}({', '.join(formatted_params)})"

def solve():
    """
    Simulates the plan to find the final state and deduces the shortest goal.
    """
    # State is represented as a set of true fluents.
    # Abbreviated names are used for clarity:
    # o1 = object_1_type_0, o2 = object_2_type_0
    # fN(...) = fluent_N(...)
    state = {'f1(o1,o2)', 'f1(o2,o1)', 'f2(o1)'}
    history = [state.copy()]

    def apply_effects(current_state, add_effects, delete_effects):
        # Effects are applied using delete-then-add semantics
        current_state.difference_update(delete_effects)
        current_state.update(add_effects)
        history.append(current_state.copy())
        return current_state

    # 1. action_3(o1, o2)
    # Pre: f1(o2,o1), f2(o1), not f0(o2,o1) -> True in initial state
    state = apply_effects(state, add_effects={'f2(o1)', 'f0(o2,o1)'}, delete_effects={'f0(o1,o2)', 'f1(o2,o1)'})

    # 2. action_2(o1, o1)
    # Pre: f2(o1) -> True
    state = apply_effects(state, add_effects={'f1(o1,o1)'}, delete_effects={'f2(o1)'})
    
    # 3. action_1(o2, o1)
    # Pre: not f2(o2), f0(o2,o1) -> True
    state = apply_effects(state, add_effects={'f2(o2)', 'f0(o2,o1)', 'f1(o2,o1)'}, delete_effects={'f2(o1)', 'f1(o1,o2)'})

    # 4. action_2(o2, o2)
    # Pre: f2(o2) -> True
    state = apply_effects(state, add_effects={'f1(o2,o2)'}, delete_effects={'f2(o2)'})

    # 5. action_1(o2, o1)
    # Pre: not f2(o2), f0(o2,o1) -> True
    state = apply_effects(state, add_effects={'f2(o2)', 'f0(o2,o1)', 'f1(o2,o1)'}, delete_effects={'f2(o1)', 'f1(o1,o2)'})
    
    # 6. action_3(o2, o2)
    # Pre: f1(o2,o2), f2(o2), not f0(o2,o2) -> True
    state = apply_effects(state, add_effects={'f2(o2)', 'f0(o2,o2)'}, delete_effects={'f0(o2,o2)', 'f1(o2,o2)'})

    # 7. action_2(o1, o2)
    # Pre: f2(o2) -> True
    state = apply_effects(state, add_effects={'f1(o1,o2)'}, delete_effects={'f2(o2)'})

    # 8. action_1(o2, o2)
    # Pre: not f2(o2), f0(o2,o2) -> True
    state = apply_effects(state, add_effects={'f2(o2)', 'f0(o2,o2)', 'f1(o2,o2)'}, delete_effects={'f2(o2)', 'f1(o2,o2)'})

    # 9. action_3(o2, o1)
    # Pre: f1(o1,o2), f2(o2), not f0(o1,o2) -> True
    state = apply_effects(state, add_effects={'f2(o2)', 'f0(o1,o2)'}, delete_effects={'f0(o2,o1)', 'f1(o1,o2)'})

    # 10. action_1(o1, o2)
    # Pre: not f2(o1), f0(o1,o2) -> True
    state = apply_effects(state, add_effects={'f2(o1)', 'f0(o1,o2)', 'f1(o1,o2)'}, delete_effects={'f2(o2)', 'f1(o2,o1)'})

    # 11. action_3(o1, o1)
    # Pre: f1(o1,o1), f2(o1), not f0(o1,o1) -> True
    state = apply_effects(state, add_effects={'f2(o1)', 'f0(o1,o1)'}, delete_effects={'f0(o1,o1)', 'f1(o1,o1)'})
    
    final_state = history[-1]
    intermediate_states = history[:-1]

    # Find the shortest goal that is only satisfied at the final step
    for size in range(1, len(final_state) + 1):
        for goal_candidate in combinations(final_state, size):
            goal_set = set(goal_candidate)
            is_truly_goal = True
            for i_state in intermediate_states:
                if goal_set.issubset(i_state):
                    is_truly_goal = False
                    break
            if is_truly_goal:
                # Found the shortest goal, format and print it.
                formatted_goal = [format_fluent(f) for f in sorted(list(goal_set))]
                print("&".join(formatted_goal))
                return

solve()