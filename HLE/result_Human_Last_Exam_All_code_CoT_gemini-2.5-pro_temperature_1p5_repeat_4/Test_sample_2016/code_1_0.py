def solve():
    """
    Simulates a plan to determine the final goal state.
    """

    # Helper functions to create fluent strings consistently
    def f0(p0, p1): return f"fluent_0({p0}, {p1})"
    def f1(p0, p1): return f"fluent_1({p0}, {p1})"
    def f2(p0): return f"fluent_2({p0})"

    # Define object names for clarity
    o1 = "object_1_type_0"
    o2 = "object_2_type_0"

    # Define the actions with their preconditions and effects
    # Preconditions are a list of (fluent, required_value) tuples, where required_value is True for presence and False for absence.
    # Effects are dictionaries with 'add' and 'del' sets of fluents.
    actions = {
        "action_0": {
            "params": 1,
            "preconds": lambda p: [(f2(p[0]), True)],
            "effects": lambda p: {"add": set(), "del": {f2(p[0])}}
        },
        "action_1": {
            "params": 2,
            "preconds": lambda p: [(f2(p[0]), False), (f0(p[0], p[1]), True)],
            "effects": lambda p: {"add": {f2(p[0]), f0(p[0], p[1]), f1(p[0], p[1])}, "del": {f2(p[1]), f1(p[1], p[0])}}
        },
        "action_2": {
            "params": 2,
            "preconds": lambda p: [(f2(p[1]), True)],
            "effects": lambda p: {"add": {f1(p[0], p[1])}, "del": {f2(p[1])}}
        },
        "action_3": {
            "params": 2,
            "preconds": lambda p: [(f1(p[1], p[0]), True), (f2(p[0]), True), (f0(p[1], p[0]), False)],
            "effects": lambda p: {"add": {f2(p[0]), f0(p[1], p[0])}, "del": {f0(p[0], p[1]), f1(p[1], p[0])}}
        }
    }

    # Set the initial state of the world
    current_state = {
        f1(o1, o2),
        f1(o2, o1),
        f2(o1)
    }

    # The optimal plan provided by the user
    plan = [
        ("action_3", (o1, o2)),
        ("action_2", (o1, o1)),
        ("action_1", (o2, o1)),
        ("action_2", (o2, o2)),
        ("action_1", (o2, o1)),
        ("action_3", (o2, o2)),
        ("action_2", (o1, o2)),
        ("action_1", (o2, o2)),
        ("action_3", (o2, o1)),
        ("action_1", (o1, o2)),
        ("action_3", (o1, o1)),
    ]

    # This set will store all fluents used to satisfy a precondition
    used_preconditions = set()

    # Execute the plan step by step
    for action_name, params in plan:
        action = actions[action_name]
        
        # Check all preconditions for the current action
        is_valid = True
        action_preconds = action["preconds"](params)
        for fluent, required_value in action_preconds:
            is_present = fluent in current_state
            if is_present != required_value:
                # This block should not be reached if the plan is valid
                print(f"Error: Precondition failure for {action_name}{params}. Required {fluent} to be {required_value}.")
                is_valid = False
                break
            # If a fluent is required to be true, we record it as an instrumental goal.
            if required_value:
                used_preconditions.add(fluent)
                
        if not is_valid:
            return # Halt execution if an error is found

        # Apply the action's effects to the state (deletions first, then additions)
        effects = action["effects"](params)
        current_state.difference_update(effects["del"])
        current_state.update(effects["add"])

    # The final state is the set of true fluents after the last action
    final_state = current_state

    # The true goal is what remains in the final state after removing instrumental goals
    goal_fluents = final_state - used_preconditions

    # Sort the fluents for a consistent, ordered output
    sorted_goal = sorted(list(goal_fluents))
    
    # Print the final goal in the specified format
    print("&".join(sorted_goal))

solve()