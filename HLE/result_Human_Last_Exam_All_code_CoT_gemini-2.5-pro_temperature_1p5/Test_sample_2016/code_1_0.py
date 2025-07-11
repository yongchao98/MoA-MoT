def solve():
    """
    Simulates a planning problem to determine the final goal state.

    The function starts with a predefined initial state and a plan (a sequence of actions).
    It iteratively applies each action to the state, verifying preconditions and
    updating the state based on the action's effects using add-after-delete semantics.
    Finally, it prints the resulting state, which is the goal.
    """

    # Initial state represented as a set of true fluents
    state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)'
    }

    # The sequence of actions in the plan
    plan = [
        ("action_3", ["object_1_type_0", "object_2_type_0"]),
        ("action_2", ["object_1_type_0", "object_1_type_0"]),
        ("action_1", ["object_2_type_0", "object_1_type_0"]),
        ("action_2", ["object_2_type_0", "object_2_type_0"]),
        ("action_1", ["object_2_type_0", "object_1_type_0"]),
        ("action_3", ["object_2_type_0", "object_2_type_0"]),
        ("action_2", ["object_1_type_0", "object_2_type_0"]),
        ("action_1", ["object_2_type_0", "object_2_type_0"]),
        ("action_3", ["object_2_type_0", "object_1_type_0"]),
        ("action_1", ["object_1_type_0", "object_2_type_0"]),
        ("action_3", ["object_1_type_0", "object_1_type_0"])
    ]

    # A helper function to format fluent names
    def f(name, params):
        return f"{name}({', '.join(params)})"

    # Process each action in the plan
    for action_name, params in plan:
        preconditions_met = False
        delete_effects = set()
        add_effects = set()

        if action_name == "action_0":
            p0 = params[0]
            if f('fluent_2', [p0]) in state:
                preconditions_met = True
                delete_effects.add(f('fluent_2', [p0]))

        elif action_name == "action_1":
            p0, p1 = params
            if f('fluent_2', [p0]) not in state and f('fluent_0', [p0, p1]) in state:
                preconditions_met = True
                add_effects.add(f('fluent_2', [p0]))
                add_effects.add(f('fluent_0', [p0, p1]))
                add_effects.add(f('fluent_1', [p0, p1]))
                delete_effects.add(f('fluent_2', [p1]))
                delete_effects.add(f('fluent_1', [p1, p0]))

        elif action_name == "action_2":
            p0, p1 = params
            if f('fluent_2', [p1]) in state:
                preconditions_met = True
                add_effects.add(f('fluent_1', [p0, p1]))
                delete_effects.add(f('fluent_2', [p1]))

        elif action_name == "action_3":
            p0, p1 = params
            if (f('fluent_1', [p1, p0]) in state and
                    f('fluent_2', [p0]) in state and
                    f('fluent_0', [p1, p0]) not in state):
                preconditions_met = True
                add_effects.add(f('fluent_2', [p0]))
                add_effects.add(f('fluent_0', [p1, p0]))
                delete_effects.add(f('fluent_0', [p0, p1]))
                delete_effects.add(f('fluent_1', [p1, p0]))

        if not preconditions_met:
            raise Exception(f"Precondition failed for {action_name}{params} in state {state}")

        # Apply effects: deletes first, then adds
        state.difference_update(delete_effects)
        state.update(add_effects)

    # Sort the final fluents for consistent output and print the goal
    goal_fluents = sorted(list(state))
    print("&".join(goal_fluents))

solve()