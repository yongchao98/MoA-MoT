import collections

class PlanningSimulator:
    def __init__(self):
        self.state = collections.defaultdict(bool)
        self.objects = ['object_1_type_0', 'object_2_type_0']
        self.initialize_state()

    def initialize_state(self):
        """Sets the initial state of the fluents."""
        self.set_fluent('fluent_1', ['object_1_type_0', 'object_2_type_0'], True)
        self.set_fluent('fluent_1', ['object_2_type_0', 'object_1_type_0'], True)
        self.set_fluent('fluent_2', ['object_1_type_0'], True)

    def get_fluent_key(self, name, params):
        """Generates a consistent key for a fluent."""
        return f"{name}({', '.join(params)})"

    def set_fluent(self, name, params, value):
        """Sets the value of a fluent in the state."""
        key = self.get_fluent_key(name, params)
        self.state[key] = value

    def get_fluent(self, name, params):
        """Gets the value of a fluent from the state."""
        key = self.get_fluent_key(name, params)
        return self.state[key]

    def action_0(self, p0):
        # Preconditions
        if not self.get_fluent('fluent_2', [p0]):
            raise ValueError(f"Precondition failed for action_0({p0}): fluent_2({p0}) is not true.")
        # Effects
        self.set_fluent('fluent_2', [p0], False)

    def action_1(self, p0, p1):
        # Preconditions
        if self.get_fluent('fluent_2', [p0]):
            raise ValueError(f"Precondition failed for action_1({p0}, {p1}): fluent_2({p0}) is true.")
        if not self.get_fluent('fluent_0', [p0, p1]):
            raise ValueError(f"Precondition failed for action_1({p0}, {p1}): fluent_0({p0}, {p1}) is not true.")
        # Effects (del then add, so add overrides)
        self.set_fluent('fluent_2', [p1], False)
        self.set_fluent('fluent_1', [p1, p0], False)
        self.set_fluent('fluent_2', [p0], True)
        self.set_fluent('fluent_0', [p0, p1], True)
        self.set_fluent('fluent_1', [p0, p1], True)

    def action_2(self, p0, p1):
        # Preconditions
        if not self.get_fluent('fluent_2', [p1]):
            raise ValueError(f"Precondition failed for action_2({p0}, {p1}): fluent_2({p1}) is not true.")
        # Effects
        self.set_fluent('fluent_2', [p1], False)
        self.set_fluent('fluent_1', [p0, p1], True)

    def action_3(self, p0, p1):
        # Preconditions
        if not self.get_fluent('fluent_1', [p1, p0]):
            raise ValueError(f"Precondition failed for action_3({p0}, {p1}): fluent_1({p1}, {p0}) is not true.")
        if not self.get_fluent('fluent_2', [p0]):
            raise ValueError(f"Precondition failed for action_3({p0}, {p1}): fluent_2({p0}) is not true.")
        if self.get_fluent('fluent_0', [p1, p0]):
            raise ValueError(f"Precondition failed for action_3({p0}, {p1}): fluent_0({p1}, {p0}) is true.")
        # Effects (del then add, so add overrides)
        self.set_fluent('fluent_0', [p0, p1], False)
        self.set_fluent('fluent_1', [p1, p0], False)
        self.set_fluent('fluent_2', [p0], True)
        self.set_fluent('fluent_0', [p1, p0], True)

    def execute_plan(self):
        """Executes the full plan."""
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
            ('action_3', ['object_1_type_0', 'object_1_type_0']),
        ]

        for action_name, params in plan:
            action_method = getattr(self, action_name)
            action_method(*params)

    def print_goal_state(self):
        """Prints the final goal state."""
        goal_fluents = [key for key, value in self.state.items() if value]
        goal_fluents.sort()
        print("&".join(goal_fluents))

if __name__ == "__main__":
    simulator = PlanningSimulator()
    simulator.execute_plan()
    simulator.print_goal_state()