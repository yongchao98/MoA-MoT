def solve():
    """
    Simulates the provided plan to determine the final goal state.
    """

    # Define the objects for easier reference
    o1 = 'object_1_type_0'
    o2 = 'object_2_type_0'

    # Initial state of the world as a set of true fluents
    current_state = {
        f'fluent_1({o1}, {o2})',
        f'fluent_1({o2}, {o1})',
        f'fluent_2({o1})'
    }

    # The plan to be executed
    plan = [
        ('action_3', o1, o2),
        ('action_2', o1, o1),
        ('action_1', o2, o1),
        ('action_2', o2, o2),
        ('action_1', o2, o1),
        ('action_3', o2, o2),
        ('action_2', o1, o2),
        ('action_1', o2, o2),
        ('action_3', o2, o1),
        ('action_1', o1, o2),
        ('action_3', o1, o1),
    ]

    # Process each action in the plan
    for action_name, p0, *params in plan:
        p1 = params[0] if params else None
        
        # A copy to modify, to avoid issues with iterating and modifying
        next_state = current_state.copy()

        if action_name == 'action_0':
            # Preconditions: fluent_2(p0)
            # Effects: not fluent_2(p0)
            if f'fluent_2({p0})' in current_state:
                next_state.remove(f'fluent_2({p0})')
        
        elif action_name == 'action_1':
            # Preconditions: not fluent_2(p0), fluent_0(p0, p1)
            # Effects add: fluent_2(p0), fluent_0(p0, p1), fluent_1(p0, p1)
            # Effects del: fluent_2(p1), fluent_1(p1, p0)
            if f'fluent_2({p0})' not in current_state and f'fluent_0({p0}, {p1})' in current_state:
                # Add effects
                next_state.add(f'fluent_2({p0})')
                next_state.add(f'fluent_0({p0}, {p1})')
                next_state.add(f'fluent_1({p0}, {p1})')
                # Del effects
                next_state.discard(f'fluent_2({p1})')
                next_state.discard(f'fluent_1({p1}, {p0})')

        elif action_name == 'action_2':
            # Preconditions: fluent_2(p1)
            # Effects add: fluent_1(p0, p1)
            # Effects del: fluent_2(p1)
            if f'fluent_2({p1})' in current_state:
                # Add effects
                next_state.add(f'fluent_1({p0}, {p1})')
                # Del effects
                next_state.remove(f'fluent_2({p1})')

        elif action_name == 'action_3':
            # Preconditions: fluent_1(p1, p0), fluent_2(p0), not fluent_0(p1, p0)
            # Effects add: fluent_2(p0), fluent_0(p1, p0)
            # Effects del: fluent_0(p0, p1), fluent_1(p1, p0)
            if (f'fluent_1({p1}, {p0})' in current_state and 
                f'fluent_2({p0})' in current_state and 
                f'fluent_0({p1}, {p0})' not in current_state):
                # Add effects
                next_state.add(f'fluent_2({p0})')
                next_state.add(f'fluent_0({p1}, {p0})')
                # Del effects
                next_state.discard(f'fluent_0({p0}, {p1})')
                next_state.discard(f'fluent_1({p1}, {p0})')

        # Update the state for the next step
        current_state = next_state
        
    # Sort the final fluents for a consistent output and print the goal
    goal_fluents = sorted(list(current_state))
    print("&".join(goal_fluents))

solve()