import collections

def solve():
    """
    Simulates a planning problem to find the goal state.
    """
    objects = ['object_1_type_0', 'object_2_type_0']
    
    # Initialize all possible fluents to False
    state = collections.defaultdict(bool)
    
    # Set initial conditions
    state['fluent_1(object_1_type_0, object_2_type_0)'] = True
    state['fluent_1(object_2_type_0, object_1_type_0)'] = True
    state['fluent_2(object_1_type_0)'] = True

    def apply_action(action_name, params, current_state):
        """
        Applies a single action to the state.
        Handles conflicting effects by letting 'true' effects overwrite 'false' effects.
        """
        p = params
        new_state = current_state.copy()

        # Define effects for each action
        effects = {
            'action_0': {'false': [f'fluent_2({p[0]})']},
            'action_1': {
                'true': [f'fluent_2({p[0]})', f'fluent_0({p[0]}, {p[1]})', f'fluent_1({p[0]}, {p[1]})'],
                'false': [f'fluent_2({p[1]})', f'fluent_1({p[1]}, {p[0]})']
            },
            'action_2': {
                'true': [f'fluent_1({p[0]}, {p[1]})'],
                'false': [f'fluent_2({p[1]})']
            },
            'action_3': {
                'true': [f'fluent_2({p[0]})', f'fluent_0({p[1]}, {p[0]})'],
                'false': [f'fluent_0({p[0]}, {p[1]})', f'fluent_1({p[1]}, {p[0]})']
            }
        }
        
        action_effects = effects.get(action_name, {})
        
        # Apply false effects first
        for fluent in action_effects.get('false', []):
            new_state[fluent] = False
            
        # Apply true effects second, potentially overwriting false ones
        for fluent in action_effects.get('true', []):
            new_state[fluent] = True
            
        return new_state

    # The provided optimal plan
    plan = [
        ('action_3', ['object_1_type_0', 'object_2_type_0']),
        ('action_2', ['object_1_type_0', 'object_1_type_0']),
        ('action_1', ['object_2_type_0', 'object_1_type_0']),
        ('action_2', ['object_2_type_0', 'object_2_type_0']),
        ('action_1', ['object_2_type_0', 'object_1_type_0']),
        ('action_3', ['object_2_type_0', 'object_2_type_0']),
        ('action_2', ['object_1_type_0', 'object_2_type_0']),
        ('action_1', ['object_2_type_0', 'object_2_type_0']),
        ('action_3', ['object_2_type_0', 'object_1_type_0']),
        ('action_1', ['object_1_type_0', 'object_2_type_0']),
        ('action_3', ['object_1_type_0', 'object_1_type_0'])
    ]

    # Execute the plan
    for action_name, params in plan:
        state = apply_action(action_name, params, state)

    # Find all true fluents in the final state
    goal_fluents = [f for f, v in sorted(state.items()) if v]
    
    # Format and print the goal
    goal_string = "&".join(goal_fluents)
    print(goal_string)

solve()