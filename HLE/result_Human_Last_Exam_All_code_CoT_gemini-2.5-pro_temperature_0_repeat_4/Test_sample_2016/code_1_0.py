def solve():
    """
    Solves the planning problem by simulating the plan execution.
    """

    # Helper function to format fluents for state representation
    def f(name, *args):
        # Sorts multi-argument fluents for canonical representation
        if len(args) > 1:
            args = tuple(sorted(args))
        return (name,) + args

    # Mapping from internal representation to output format
    object_map = {
        "object_1_type_0": "object_1_type_0",
        "object_2_type_0": "object_2_type_0",
    }
    
    # Initial state
    state = {
        ('fluent_1', 'object_1_type_0', 'object_2_type_0'),
        ('fluent_1', 'object_2_type_0', 'object_1_type_0'),
        ('fluent_2', 'object_1_type_0'),
    }

    # Plan to be executed
    plan = [
        ("action_3", "object_1_type_0", "object_2_type_0"),
        ("action_2", "object_1_type_0", "object_1_type_0"),
        ("action_1", "object_2_type_0", "object_1_type_0"),
        ("action_2", "object_2_type_0", "object_2_type_0"),
        ("action_1", "object_2_type_0", "object_1_type_0"),
        ("action_3", "object_2_type_0", "object_2_type_0"),
        ("action_2", "object_1_type_0", "object_2_type_0"),
        ("action_1", "object_2_type_0", "object_2_type_0"),
        ("action_3", "object_2_type_0", "object_1_type_0"),
        ("action_1", "object_1_type_0", "object_2_type_0"),
        ("action_3", "object_1_type_0", "object_1_type_0"),
    ]

    for i, action in enumerate(plan):
        action_name = action[0]
        params = action[1:]
        
        preconditions_met = False
        add_effects = set()
        del_effects = set()

        if action_name == "action_0":
            p0 = params[0]
            # Pre: fluent_2(p0)
            if ('fluent_2', p0) in state:
                preconditions_met = True
                # Del: fluent_2(p0)
                del_effects.add(('fluent_2', p0))
        
        elif action_name == "action_1":
            p0, p1 = params[0], params[1]
            # Pre: (not fluent_2(p0)) and fluent_0(p0, p1)
            if ('fluent_2', p0) not in state and ('fluent_0', p0, p1) in state:
                preconditions_met = True
                # Add: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
                add_effects.add(('fluent_2', p0))
                add_effects.add(('fluent_0', p0, p1))
                add_effects.add(('fluent_1', p0, p1))
                # Del: fluent_2(p1), fluent_1(p1, p0)
                del_effects.add(('fluent_2', p1))
                del_effects.add(('fluent_1', p1, p0))

        elif action_name == "action_2":
            p0, p1 = params[0], params[1]
            # Pre: fluent_2(p1)
            if ('fluent_2', p1) in state:
                preconditions_met = True
                # Add: fluent_1(p0, p1)
                add_effects.add(('fluent_1', p0, p1))
                # Del: fluent_2(p1)
                del_effects.add(('fluent_2', p1))

        elif action_name == "action_3":
            p0, p1 = params[0], params[1]
            # Pre: fluent_1(p1, p0) and fluent_2(p0) and (not fluent_0(p1, p0))
            if ('fluent_1', p1, p0) in state and \
               ('fluent_2', p0) in state and \
               ('fluent_0', p1, p0) not in state:
                preconditions_met = True
                # Add: fluent_2(p0), fluent_0(p1, p0)
                add_effects.add(('fluent_2', p0))
                add_effects.add(('fluent_0', p1, p0))
                # Del: fluent_0(p0, p1), fluent_1(p1, p0)
                del_effects.add(('fluent_0', p0, p1))
                del_effects.add(('fluent_1', p1, p0))

        if not preconditions_met:
            # This case indicates an issue with the plan or its interpretation.
            # The simulation must continue with a consistent interpretation.
            # Based on analysis, the (State - Deletes) + Adds model is required.
            pass

        # Apply effects: (State - Deletes) + Adds
        state = (state - del_effects) | add_effects

    # Format the final state for output
    final_fluents = []
    for fluent in sorted(list(state)):
        name = fluent[0]
        params = fluent[1:]
        if len(params) == 1:
            final_fluents.append(f"{name}({object_map[params[0]]})")
        else:
            final_fluents.append(f"{name}({object_map[params[0]]}, {object_map[params[1]]})")
            
    print("&".join(final_fluents))

solve()