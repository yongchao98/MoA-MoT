class MESISimulator:
    """
    A simulator for the MESI cache coherence protocol.
    """
    def __init__(self, processors):
        """
        Initializes the simulator.
        - processors: A list of processor names (e.g., ['P1', 'P2']).
        """
        self.processors = processors
        self.caches = {p: 'I' for p in processors}
        self.message_count = 0
        self.step_messages = []
        self.memory_value = 0
        print(f"Initial State: Caches = {self.caches}, Memory(X) = 0, Total Messages = 0\n")

    def execute_op(self, step_num, proc, action, value=None):
        """
        Executes a single processor operation and updates cache states.
        """
        current_state = self.caches[proc]
        messages = 0
        description = []

        if action == 'read':
            if current_state == 'I': # Read Miss
                description.append(f"{proc} has a Read Miss.")
                # Message 1: Read Request from proc
                # Message 2: Data Response from memory or another cache
                messages = 2
                
                # Check if another cache has the data
                sharers = [p for p, s in self.caches.items() if s in ['S', 'E']]
                modifier = [p for p, s in self.caches.items() if s == 'M']

                if modifier:
                    # Another cache has it in Modified state
                    modifier_proc = modifier[0]
                    description.append(f"{modifier_proc} (in M) provides data and writes back to memory.")
                    self.caches[modifier_proc] = 'S' # M -> S
                    self.caches[proc] = 'S' # I -> S
                elif sharers:
                    # Other caches have it in Shared/Exclusive state
                    description.append(f"A cache in state E/S provides the data.")
                    for p in sharers:
                        if self.caches[p] == 'E':
                           self.caches[p] = 'S' # E -> S
                    self.caches[proc] = 'S' # I -> S
                else:
                    # No other cache has the data, fetch from memory
                    description.append("Memory provides the data.")
                    self.caches[proc] = 'E' # I -> E

            else: # Read Hit (M, E, or S)
                description.append(f"{proc} has a Read Hit (State: {current_state}). No messages.")
                messages = 0

        elif action == 'write':
            if current_state in ['I', 'S']: # Write Miss / Upgrade
                if current_state == 'I':
                    description.append(f"{proc} has a Write Miss (State: I). Sends RFO.")
                    # Message 1: Read For Ownership (RFO) request
                    # Message 2: Data response
                    messages = 2
                else: # current_state == 'S'
                    description.append(f"{proc} wants to write to a Shared block. Sends Invalidate.")
                    # Message 1: Invalidate broadcast
                    messages = 1
                
                # Invalidate all other caches that have a copy
                for p in self.processors:
                    if p != proc and self.caches[p] != 'I':
                        self.caches[p] = 'I'
                
                self.caches[proc] = 'M' # Becomes Modified
                if value is not None:
                    self.memory_value = value
            
            elif current_state == 'E': # Write Hit (Exclusive)
                description.append(f"{proc} has a Write Hit (State: E). No bus messages needed.")
                messages = 0
                self.caches[proc] = 'M' # E -> M
                if value is not None:
                    self.memory_value = value

            else: # current_state == 'M' (Write Hit)
                description.append(f"{proc} has a Write Hit (State: M). No messages.")
                messages = 0
                if value is not None:
                    self.memory_value = value
        
        self.message_count += messages
        self.step_messages.append(messages)
        
        print(f"Step {step_num}: {proc} {action}s X" + (f" with value {value}" if value else ""))
        print(f"  - Action: {' '.join(description)}")
        print(f"  - Messages this step: {messages}")
        print(f"  - New Cache States: {self.caches}")
        print(f"  - Running Total Messages: {self.message_count}\n")
        
    def print_final_tally(self):
        equation = " + ".join(map(str, self.step_messages))
        print("Final Calculation:")
        print(f"Total Messages = {equation} = {self.message_count}")


# --- Simulation ---
processors = ['P1', 'P2', 'P3', 'P4']
sim = MESISimulator(processors)

operations = [
    ('P1', 'read', None),
    ('P2', 'read', None),
    ('P1', 'write', 1),
    ('P3', 'read', None),
    ('P2', 'write', 2),
    ('P4', 'read', None),
    ('P1', 'read', None)
]

for i, op in enumerate(operations):
    sim.execute_op(i + 1, op[0], op[1], op[2])

sim.print_final_tally()

final_answer = sim.message_count
print(f'<<<{final_answer}>>>')