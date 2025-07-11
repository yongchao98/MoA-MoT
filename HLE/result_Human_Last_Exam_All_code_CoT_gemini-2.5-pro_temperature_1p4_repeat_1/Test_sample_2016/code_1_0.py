import collections

class State:
    """
    A class to represent the state of the planning problem.
    It holds the set of true and false fluents.
    """
    def __init__(self, objects):
        self.fluents = {}
        # Enumerate all possible fluents and set them to False initially
        for p1 in objects:
            self.fluents[f'fluent_2({p1})'] = False
            for p2 in objects:
                self.fluents[f'fluent_0({p1}, {p2})'] = False
                self.fluents[f'fluent_1({p1}, {p2})'] = False

    def set_initial_conditions(self, initial_fluents):
        """Sets the initial true fluents."""
        for fluent in initial_fluents:
            self.fluents[fluent] = True

    def check_preconditions(self, preconditions):
        """Checks if all preconditions are met in the current state."""
        for fluent, value in preconditions.items():
            if self.fluents.get(fluent) != value:
                # Uncomment the following line for debugging precondition failures
                # print(f"Precondition failed: {fluent} is not {value}")
                return False
        return True

    def apply_effects(self, effects):
        """Applies the effects of an action to the state (add-after-delete)."""
        new_state = self.fluents.copy()
        
        # Apply delete effects first
        for fluent in effects.get('delete', []):
            new_state[fluent] = False
            
        # Apply add effects
        for fluent in effects.get('add', []):
            new_state[fluent] = True
            
        self.fluents = new_state

    def get_true_fluents(self):
        """Returns a sorted list of fluents that are currently true."""
        return sorted([k for k, v in self.fluents.items() if v])

# --- Action Definitions ---

def action_0(state, p0):
    preconditions = {f'fluent_2({p0})': True}
    if not state.check_preconditions(preconditions):
        raise ValueError(f"Preconditions not met for action_0({p0})")
    
    effects = {'delete': [f'fluent_2({p0})']}
    state.apply_effects(effects)

def action_1(state, p0, p1):
    preconditions = {
        f'fluent_2({p0})': False,
        f'fluent_0({p0}, {p1})': True
    }
    if not state.check_preconditions(preconditions):
        raise ValueError(f"Preconditions not met for action_1({p0}, {p1})")
        
    effects = {
        'add': [f'fluent_2({p0})', f'fluent_0({p0}, {p1})', f'fluent_1({p0}, {p1})'],
        'delete': [f'fluent_2({p1})', f'fluent_1({p1}, {p0})']
    }
    state.apply_effects(effects)

def action_2(state, p0, p1):
    preconditions = {f'fluent_2({p1})': True}
    if not state.check_preconditions(preconditions):
        raise ValueError(f"Preconditions not met for action_2({p0}, {p1})")
        
    effects = {
        'add': [f'fluent_1({p0}, {p1})'],
        'delete': [f'fluent_2({p1})']
    }
    state.apply_effects(effects)

def action_3(state, p0, p1):
    preconditions = {
        f'fluent_1({p1}, {p0})': True,
        f'fluent_2({p0})': True,
        f'fluent_0({p1}, {p0})': False
    }
    if not state.check_preconditions(preconditions):
        raise ValueError(f"Preconditions not met for action_3({p0}, {p1})")

    effects = {
        'add': [f'fluent_2({p0})', f'fluent_0({p1}, {p0})'],
        'delete': [f'fluent_0({p0}, {p1})', f'fluent_1({p1}, {p0})']
    }
    state.apply_effects(effects)

def solve():
    """
    Main function to simulate the plan and find the goal.
    """
    # Define objects with their full names from the problem
    obj_map = {
        'obj1': 'object_1_type_0',
        'obj2': 'object_2_type_0'
    }
    objects = list(obj_map.values())
    
    # Initialize state
    world_state = State(objects)
    
    # Set initial conditions
    initial_conditions = [
        f'fluent_1({obj_map["obj1"]}, {obj_map["obj2"]})',
        f'fluent_1({obj_map["obj2"]}, {obj_map["obj1"]})',
        f'fluent_2({obj_map["obj1"]})'
    ]
    world_state.set_initial_conditions(initial_conditions)

    # The optimal plan
    plan = [
        lambda s: action_3(s, obj_map['obj1'], obj_map['obj2']),
        lambda s: action_2(s, obj_map['obj1'], obj_map['obj1']),
        lambda s: action_1(s, obj_map['obj2'], obj_map['obj1']),
        lambda s: action_2(s, obj_map['obj2'], obj_map['obj2']),
        lambda s: action_1(s, obj_map['obj2'], obj_map['obj1']),
        lambda s: action_3(s, obj_map['obj2'], obj_map['obj2']),
        lambda s: action_2(s, obj_map['obj1'], obj_map['obj2']),
        lambda s: action_1(s, obj_map['obj2'], obj_map['obj2']),
        lambda s: action_3(s, obj_map['obj2'], obj_map['obj1']),
        lambda s: action_1(s, obj_map['obj1'], obj_map['obj2']),
        lambda s: action_3(s, obj_map['obj1'], obj_map['obj1']),
    ]

    # Execute the plan
    for action in plan:
        action(world_state)
        
    # Get the final state and format as the goal
    goal_fluents = world_state.get_true_fluents()
    goal_string = "&".join(goal_fluents)
    
    print(goal_string)

solve()