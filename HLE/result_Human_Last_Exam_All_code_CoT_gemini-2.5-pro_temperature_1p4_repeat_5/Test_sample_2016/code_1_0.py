import collections

def solve():
    """
    Simulates a planning problem to determine the goal state.
    """

    # Step 1 & 2: Represent and Initialize the World State
    objects = {
        'o1': 'object_1_type_0',
        'o2': 'object_2_type_0'
    }

    state = {
        # fluent_0
        'f0_o1_o1': False, 'f0_o1_o2': False, 'f0_o2_o1': False, 'f0_o2_o2': False,
        # fluent_1
        'f1_o1_o1': False, 'f1_o1_o2': True,  'f1_o2_o1': True,  'f1_o2_o2': False,
        # fluent_2
        'f2_o1': True, 'f2_o2': False
    }

    # Step 3: Implement the Actions
    
    def apply_action(action_name, params, current_state):
        """
        Applies a given action to the current state.
        This function implements an "add wins" semantic for conflicting effects.
        """
        p_len = len(params)
        new_state = current_state.copy()
        
        preconditions_met = False
        add_effects = []
        del_effects = []

        if action_name == 'action_0':
            p0, = params
            if current_state[f'f2_{p0}']:
                preconditions_met = True
                del_effects.append(f'f2_{p0}')

        elif action_name == 'action_1':
            p0, p1 = params
            if not current_state[f'f2_{p0}'] and current_state[f'f0_{p0}_{p1}']:
                preconditions_met = True
                add_effects.extend([f'f2_{p0}', f'f0_{p0}_{p1}', f'f1_{p0}_{p1}'])
                del_effects.extend([f'f2_{p1}', f'f1_{p1}_{p0}'])

        elif action_name == 'action_2':
            p0, p1 = params
            if current_state[f'f2_{p1}']:
                preconditions_met = True
                add_effects.append(f'f1_{p0}_{p1}')
                del_effects.append(f'f2_{p1}')

        elif action_name == 'action_3':
            p0, p1 = params
            if current_state[f'f1_{p1}_{p0}'] and current_state[f'f2_{p0}'] and not current_state[f'f0_{p1}_{p0}']:
                preconditions_met = True
                add_effects.extend([f'f2_{p0}', f'f0_{p1}_{p0}'])
                del_effects.extend([f'f0_{p0}_{p1}', f'f1_{p1}_{p0}'])
        
        if not preconditions_met:
            raise ValueError(f"Preconditions for {action_name}{params} not met.")

        # Apply effects: "add wins" by applying deletes first, then adds.
        for effect in del_effects:
            new_state[effect] = False
        for effect in add_effects:
            new_state[effect] = True
            
        return new_state

    # Step 4: Simulate the Plan
    plan = [
        ('action_3', ('o1', 'o2')),
        ('action_2', ('o1', 'o1')),
        ('action_1', ('o2', 'o1')),
        ('action_2', ('o2', 'o2')),
        ('action_1', ('o2', 'o1')),
        ('action_3', ('o2', 'o2')),
        ('action_2', ('o1', 'o2')),
        ('action_1', ('o2', 'o2')),
        ('action_3', ('o2', 'o1')),
        ('action_1', ('o1', 'o2')),
        ('action_3', ('o1', 'o1')),
    ]

    for i, (action, params) in enumerate(plan):
        state = apply_action(action, params, state)

    # Step 5 & 6: Identify and Format the Goal
    final_true_fluents = []
    # Sort keys for consistent output order
    for fluent_key in sorted(state.keys()):
        if state[fluent_key]:
            parts = fluent_key.split('_')
            fluent_name = f"fluent_{parts[0][1:]}"
            if len(parts) == 2:
                param1 = objects[parts[1]]
                final_true_fluents.append(f"{fluent_name}({param1})")
            elif len(parts) == 3:
                param1 = objects[parts[1]]
                param2 = objects[parts[2]]
                final_true_fluents.append(f"{fluent_name}({param1}, {param2})")
    
    print("&".join(final_true_fluents))

solve()