class MESISimulator:
    """
    Simulates the MESI protocol for a sequence of operations.
    Counts the cache coherence messages exchanged. A message is a bus
    transaction like a request or a data response.
    - Read/Write Miss (BusRd, BusRdX): 2 messages (request + data response)
    - Write Upgrade (BusUpgr): 1 message (invalidate broadcast)
    """

    def __init__(self, processors):
        """Initializes the simulator."""
        self.processors = processors
        self.caches = {p: {'state': 'I', 'value': None} for p in self.processors}
        self.memory_value = 0
        self.total_messages = 0
        self.step_count = 0
        print(f"Initial State: X = 0 in memory. All caches are in 'Invalid' state.\n")

    def _get_sharers(self, requester):
        """Finds other processors sharing the cache line."""
        return [p for p, cache in self.caches.items() if p != requester and cache['state'] == 'S']

    def _get_exclusive_or_modified_holder(self, requester):
        """Finds if any other processor has the line in E or M state."""
        for p, cache in self.caches.items():
            if p != requester and cache['state'] in ['E', 'M']:
                return p
        return None

    def read(self, processor):
        """Simulates a read operation by a processor."""
        self.step_count += 1
        op_string = f"{processor} reads X"
        current_cache = self.caches[processor]
        messages_this_step = 0
        explanation = []

        if current_cache['state'] in ['M', 'E', 'S']:
            # Read Hit: No messages, no state changes
            explanation.append(f"{processor} has a read hit. Data is read locally.")
        else: # State is 'I' -> Read Miss
            requester_state = 'I'
            messages_this_step = 2  # BusRd request + Data response
            
            # Check other caches
            other_valid_holder = self._get_exclusive_or_modified_holder(processor)
            sharers = self._get_sharers(processor)

            if other_valid_holder:
                holder_cache = self.caches[other_valid_holder]
                explanation.append(f"{processor} has a read miss and sends a BusRd request.")
                explanation.append(f"{other_valid_holder} (in state '{holder_cache['state']}') snoops the request and provides the data.")
                explanation.append("This is a cache-to-cache transfer.")
                
                # If state was 'M', data is also written back to memory
                if holder_cache['state'] == 'M':
                    explanation.append(f"{other_valid_holder} also writes its modified data back to memory.")
                    self.memory_value = holder_cache['value']
                
                # Both caches become Shared
                holder_cache['state'] = 'S'
                current_cache['state'] = 'S'
                current_cache['value'] = holder_cache['value']
            elif sharers:
                # One of the sharers provides the data
                data_supplier = sharers[0]
                supplier_cache = self.caches[data_supplier]
                explanation.append(f"{processor} has a read miss and sends a BusRd request.")
                explanation.append(f"A sharing cache ({data_supplier}) provides the data.")
                current_cache['state'] = 'S'
                current_cache['value'] = supplier_cache['value']
            else:
                # No other cache has the data, fetch from memory
                explanation.append(f"{processor} has a read miss and sends a BusRd request.")
                explanation.append("No other cache has the data, so memory responds.")
                current_cache['state'] = 'E' # Becomes Exclusive
                current_cache['value'] = self.memory_value

        self._print_step_result(op_string, messages_this_step, explanation)


    def write(self, processor, value):
        """Simulates a write operation by a processor."""
        self.step_count += 1
        op_string = f"{processor} writes X = {value}"
        current_cache = self.caches[processor]
        messages_this_step = 0
        explanation = []

        if current_cache['state'] == 'M':
            # Write Hit, no messages
            explanation.append(f"{processor} is already in 'Modified' state. Write happens locally.")
        elif current_cache['state'] == 'E':
            # Write Hit, no bus messages. Just change state to M.
            explanation.append(f"{processor} is in 'Exclusive' state. It writes locally and transitions to 'Modified'.")
            current_cache['state'] = 'M'
        elif current_cache['state'] == 'S':
            # Write Hit, but needs to invalidate others
            messages_this_step = 1 # BusUpgr (Invalidate)
            explanation.append(f"{processor} is in 'Shared' state. It sends an Invalidate message on the bus.")
            sharers = self._get_sharers(processor)
            for p in sharers:
                self.caches[p]['state'] = 'I'
                explanation.append(f"{p}'s copy is invalidated.")
            current_cache['state'] = 'M'
        else: # State is 'I' -> Write Miss
            messages_this_step = 2 # BusRdX request + data response
            explanation.append(f"{processor} has a write miss and sends a BusRdX (Read with Intent to Modify) request.")
            
            # Invalidate all other copies and get data
            for p, cache in self.caches.items():
                if p != processor and cache['state'] != 'I':
                    explanation.append(f"{p} (in state '{cache['state']}') snoops the request, provides data if required, and invalidates its copy.")
                    if cache['state'] == 'M':
                        self.memory_value = cache['value']
                    cache['state'] = 'I'
            
            current_cache['state'] = 'M'

        current_cache['value'] = value
        self._print_step_result(op_string, messages_this_step, explanation)

    def _print_step_result(self, op_string, messages, explanation):
        """Prints the results of a single simulation step."""
        print(f"--- Step {self.step_count}: {op_string} ---")
        for line in explanation:
            print(f"  - {line}")
        
        self.total_messages += messages
        print(f"\n  Messages this step: {messages}")
        print(f"  Cache States: " + ", ".join([f"{p}:{c['state']}" for p, c in self.caches.items()]))
        print(f"  Total Messages so far: {self.total_messages}\n")


# Main execution
if __name__ == "__main__":
    processors = ['P1', 'P2', 'P3', 'P4']
    sim = MESISimulator(processors)

    # Sequence of operations
    sim.read('P1')
    sim.read('P2')
    sim.write('P1', 1)
    sim.read('P3')
    sim.write('P2', 2)
    sim.read('P4')
    sim.read('P1')

    print("=" * 30)
    print("Final Calculation:")
    print("The final total number of cache coherence messages is the sum of messages from each step.")
    print("Equation: 2 + 2 + 1 + 2 + 2 + 2 + 2")
    print(f"Total messages = {sim.total_messages}")
    print("=" * 30)
    print(f"<<<{sim.total_messages}>>>")