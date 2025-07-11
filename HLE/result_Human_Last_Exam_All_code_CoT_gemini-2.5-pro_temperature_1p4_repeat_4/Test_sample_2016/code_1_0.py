def solve():
    """
    Simulates a plan to find the final goal state.
    """
    # Define the objects
    objects = ["object_1_type_0", "object_2_type_0"]
    o1, o2 = objects[0], objects[1]

    # --- State Representation ---
    # A single dictionary to hold the state of all fluents.
    # Keys are tuples: (fluent_name, param1, ...)
    # Values are booleans.
    state = {}

    # --- State Initialization ---
    # Initialize all possible fluents to False
    for p1 in objects:
        state[('fluent_2', p1)] = False
        for p2 in objects:
            state[('fluent_0', p1, p2)] = False
            state[('fluent_1', p1, p2)] = False

    # Set initial true conditions
    state[('fluent_1', o1, o2)] = True
    state[('fluent_1', o2, o1)] = True
    state[('fluent_2', o1)] = True

    # --- Action Definitions ---
    # Each action checks preconditions (for correctness) and applies effects.

    def action_0(p0):
        # Pre: fluent_2(p0)
        # Eff: not fluent_2(p0)
        if not state[('fluent_2', p0)]:
            raise ValueError(f"Precondition failed for action_0({p0})")
        state[('fluent_2', p0)] = False

    def action_1(p0, p1):
        # Pre: (not fluent_2(p0)) and fluent_0(p0, p1)
        # Add: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
        # Del: fluent_2(p1), fluent_1(p1, p0)
        if not (not state[('fluent_2', p0)] and state[('fluent_0', p0, p1)]):
            raise ValueError(f"Precondition failed for action_1({p0}, {p1})")
        # Apply del effects first
        state[('fluent_2', p1)] = False
        state[('fluent_1', p1, p0)] = False
        # Apply add effects
        state[('fluent_2', p0)] = True
        state[('fluent_0', p0, p1)] = True
        state[('fluent_1', p0, p1)] = True

    def action_2(p0, p1):
        # Pre: fluent_2(p1)
        # Add: fluent_1(p0, p1)
        # Del: fluent_2(p1)
        if not state[('fluent_2', p1)]:
            raise ValueError(f"Precondition failed for action_2({p0}, {p1})")
        state[('fluent_2', p1)] = False
        state[('fluent_1', p0, p1)] = True

    def action_3(p0, p1):
        # Pre: fluent_1(p1, p0) and fluent_2(p0) and (not fluent_0(p1, p0))
        # Add: fluent_2(p0), fluent_0(p1, p0)
        # Del: fluent_0(p0, p1), fluent_1(p1, p0)
        if not (state[('fluent_1', p1, p0)] and state[('fluent_2', p0)] and not state[('fluent_0', p1, p0)]):
            raise ValueError(f"Precondition failed for action_3({p0}, {p1})")
        # Apply del effects first
        state[('fluent_0', p0, p1)] = False
        state[('fluent_1', p1, p0)] = False
        # Apply add effects
        state[('fluent_2', p0)] = True
        state[('fluent_0', p1, p0)] = True

    # --- Plan Execution ---
    plan = [
        (action_3, (o1, o2)),
        (action_2, (o1, o1)),
        (action_1, (o2, o1)),
        (action_2, (o2, o2)),
        (action_1, (o2, o1)),
        (action_3, (o2, o2)),
        (action_2, (o1, o2)),
        (action_1, (o2, o2)),
        (action_3, (o2, o1)),
        (action_1, (o1, o2)),
        (action_3, (o1, o1)),
    ]

    for action, params in plan:
        action(*params)

    # --- Determine and Print Final State ---
    final_true_fluents = []
    # Sort keys for a consistent, deterministic output order
    for key in sorted(state.keys()):
        if state[key]:
            fluent_name = key[0]
            fluent_params = ", ".join(key[1:])
            final_true_fluents.append(f"{fluent_name}({fluent_params})")
    
    goal_string = "&".join(final_true_fluents)
    print(goal_string)

solve()
<<<fluent_0(object_1_type_0, object_1_type_0)&fluent_0(object_1_type_0, object_2_type_0)&fluent_0(object_2_type_0, object_2_type_0)&fluent_1(object_1_type_0, object_2_type_0)&fluent_1(object_2_type_0, object_2_type_0)&fluent_2(object_1_type_0)>>>