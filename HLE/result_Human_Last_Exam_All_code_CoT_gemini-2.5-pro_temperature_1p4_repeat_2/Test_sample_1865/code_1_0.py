import collections

class MesiSimulator:
    """
    Simulates the MESI protocol for a multiprocessor system.
    """
    def __init__(self, processors, initial_memory_value=0):
        self.processors = processors
        self.states = {p: 'I' for p in self.processors}
        self.memory_value = initial_memory_value
        self.message_count = 0
        self.message_log = []
        self.message_equation_parts = []
        print(f"Initial State: All caches are Invalid. Memory X = {self.memory_value}. Messages = 0.")
        print("-" * 60)

    def _print_state(self):
        """Prints the current state of all processor caches."""
        state_str = ", ".join([f"{p}: {s}" for p, s in self.states.items()])
        print(f"-> Cache States: [{state_str}]")
        print(f"-> Memory X = {self.memory_value}")

    def read(self, processor):
        """Simulates a read operation by a processor."""
        print(f"\nProcessing: {processor} reads X")
        messages_this_step = 0
        
        # Cache hit
        if self.states[processor] in ['M', 'E', 'S']:
            print(f"- {processor}'s cache has X in '{self.states[processor]}' state. Cache hit.")
            print("- No message needed.")
        
        # Cache miss
        elif self.states[processor] == 'I':
            print(f"- {processor}'s cache is Invalid. This is a Read Miss.")
            print(f"- {processor} sends a 'Read Miss' message.")
            self.message_count += 1
            messages_this_step += 1
            
            # Check other caches
            other_states = {p: s for p, s in self.states.items() if p != processor}
            
            # If another cache has it in Modified state
            if 'M' in other_states.values():
                owner = [p for p, s in other_states.items() if s == 'M'][0]
                print(f"- {owner} has the block in 'M' state. It supplies the data to {processor} and writes back to memory.")
                # The owner's state changes to Shared
                self.states[owner] = 'S'
                # The memory is updated
                # The value would be supplied by the owner cache, we don't simulate exact values here
                self.memory_value = f"Value from {owner}"
                # The requester's state becomes Shared
                self.states[processor] = 'S'
            # If another cache has it in Exclusive state
            elif 'E' in other_states.values():
                owner = [p for p, s in other_states.items() if s == 'E'][0]
                print(f"- {owner} has the block in 'E' state. It supplies the data.")
                # The owner's state changes to Shared
                self.states[owner] = 'S'
                # The requester's state becomes Shared
                self.states[processor] = 'S'
            # If other caches have it in Shared state
            elif 'S' in other_states.values():
                print("- Another cache supplies the data.")
                # The requester's state becomes Shared
                self.states[processor] = 'S'
            # If no other cache has the data
            else:
                print("- Memory supplies the data. No other cache shares the block.")
                # The requester's state becomes Exclusive
                self.states[processor] = 'E'
        
        self.message_equation_parts.append(str(messages_this_step))
        self._print_state()
        print("-" * 60)

    def write(self, processor, value):
        """Simulates a write operation by a processor."""
        print(f"\nProcessing: {processor} writes X = {value}")
        messages_this_step = 0

        # Case 1: Hit in M or E state
        if self.states[processor] in ['M', 'E']:
            print(f"- {processor}'s cache has X in '{self.states[processor]}' state. Write hit.")
            print("- No coherence message needed.")
            self.states[processor] = 'M' # State becomes/remains Modified
        
        # Case 2: Hit in S state (Upgrade needed)
        elif self.states[processor] == 'S':
            print("- {processor}'s cache has X in 'S' state. This requires an upgrade.")
            print(f"- {processor} sends an 'Invalidate' message to invalidate other shared copies.")
            self.message_count += 1
            messages_this_step += 1
            self.states[processor] = 'M'
            for p in self.processors:
                if p != processor and self.states[p] == 'S':
                    self.states[p] = 'I'
                    print(f"- {p}'s copy is invalidated.")
            self.memory_value = f"{self.memory_value} (Stale)"

        # Case 3: Miss in I state (Write Miss)
        elif self.states[processor] == 'I':
            print("- {processor}'s cache is Invalid. This is a Write Miss.")
            print("- {processor} sends a 'Read For Ownership' (RFO) message.")
            self.message_count += 1
            messages_this_step += 1
            # Invalidate all other copies
            for p in self.processors:
                if p != processor:
                    if self.states[p] in ['M', 'S', 'E']:
                        # If a cache had 'M', it would write back to memory before invalidating.
                        if self.states[p] == 'M':
                            print(f"- {p} has X in 'M' state, it writes back to memory and invalidates its copy.")
                            self.memory_value = f"Value from {p}"
                        else:
                            print(f"- {p}'s copy in state '{self.states[p]}' is invalidated.")
                        self.states[p] = 'I'

            self.states[processor] = 'M'
            self.memory_value = f"{self.memory_value} (Stale)"

        self.message_equation_parts.append(str(messages_this_step))
        self._print_state()
        print("-" * 60)
        
    def get_final_count(self):
        """Prints the final result."""
        operation_names = ["P1 Read", "P2 Read", "P1 Write", "P3 Read", "P2 Write", "P4 Read", "P1 Read"]
        equation_parts_with_ops = [f"{self.message_equation_parts[i]} ({op})" for i, op in enumerate(operation_names)]
        equation = " + ".join(equation_parts_with_ops)
        print(f"\nFinal Calculation:")
        print(f"Total Messages = {equation} = {self.message_count}")
        

# --- Main Execution ---
processors_list = ['P1', 'P2', 'P3', 'P4']
simulator = MesiSimulator(processors_list)

# Sequence of operations
simulator.read('P1')         # 1. P1 reads X
simulator.read('P2')         # 2. P2 reads X
simulator.write('P1', 1)     # 3. P1 writes X = 1
simulator.read('P3')         # 4. P3 reads X
simulator.write('P2', 2)     # 5. P2 writes X = 2
simulator.read('P4')         # 6. P4 reads X
simulator.read('P1')         # 7. P1 reads X

simulator.get_final_count()
print(f"<<<{simulator.message_count}>>>")