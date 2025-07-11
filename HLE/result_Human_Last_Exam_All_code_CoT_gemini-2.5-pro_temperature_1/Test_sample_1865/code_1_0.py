import sys

class MESISimulator:
    """
    A simulator for the MESI cache coherence protocol.
    It tracks the state of a shared variable in multiple processor caches
    and counts the coherence messages exchanged.
    """

    def __init__(self, processor_ids):
        """
        Initializes the system.
        - processor_ids: A list of processor names (e.g., ['P1', 'P2']).
        """
        self.processors = processor_ids
        # All caches start in Invalid state
        self.cache_states = {p_id: 'I' for p_id in self.processors}
        self.memory_value = 0
        self.message_count = 0
        self.message_log = []
        print("--- Initial State ---")
        self._print_status("System initialized. X=0 in memory. All caches are Invalid.")
        print("-" * 25)

    def _print_status(self, description):
        """Prints the current status of the caches and message count."""
        states_str = ", ".join([f"{p}: {s}" for p, s in self.cache_states.items()])
        print(f"Description: {description}")
        print(f"Cache States: {{ {states_str} }}")
        print(f"Cumulative Messages: {self.message_count}\n")

    def read(self, p_id):
        """Simulates a read operation by a processor."""
        operation_str = f"Step: {len(self.message_log) + 1}. {p_id} reads X"
        print(operation_str)
        
        state = self.cache_states[p_id]

        if state in ['M', 'E', 'S']:
            # Read Hit: No message needed
            self._print_status(f"{p_id} has a Read Hit. No coherence message is sent.")
        elif state == 'I':
            # Read Miss: A message is sent
            self.message_count += 1
            self.message_log.append(1)
            
            # Check states of other caches to determine action
            other_states = [s for p, s in self.cache_states.items() if p != p_id]
            
            if 'M' in other_states:
                # Another cache has the data in Modified state
                provider_p = [p for p, s in self.cache_states.items() if s == 'M'][0]
                self.cache_states[provider_p] = 'S' # Provider changes to Shared
                self.cache_states[p_id] = 'S'       # Requester becomes Shared
                desc = (f"{p_id} has a Read Miss. Sends 'Read Miss' message (1). "
                        f"{provider_p} (in M) provides data, writes back to memory, and changes to S. "
                        f"{p_id} becomes S.")
            elif 'E' in other_states:
                # Another cache has the data in Exclusive state
                provider_p = [p for p, s in self.cache_states.items() if s == 'E'][0]
                self.cache_states[provider_p] = 'S' # Provider changes to Shared
                self.cache_states[p_id] = 'S'       # Requester becomes Shared
                desc = (f"{p_id} has a Read Miss. Sends 'Read Miss' message (1). "
                        f"{provider_p} (in E) provides data and changes to S. {p_id} becomes S.")
            else:
                # Data comes from memory, or another Shared cache
                # If others are 'S', one of them provides data. If all 'I', memory provides.
                # In either case, the new state is Shared if others have it, else Exclusive.
                if 'S' in other_states:
                    self.cache_states[p_id] = 'S'
                    desc = (f"{p_id} has a Read Miss. Sends 'Read Miss' message (1). "
                            f"Another cache (in S) provides data. {p_id} becomes S.")
                else: # All others must be 'I'
                    self.cache_states[p_id] = 'E'
                    desc = (f"{p_id} has a Read Miss. Sends 'Read Miss' message (1). "
                            f"Memory provides data. {p_id} becomes Exclusive (E).")
            self._print_status(desc)
        print("-" * 25)

    def write(self, p_id, value):
        """Simulates a write operation by a processor."""
        operation_str = f"Step: {len(self.message_log) + 1}. {p_id} writes X = {value}"
        print(operation_str)
        
        state = self.cache_states[p_id]
        
        if state == 'M':
            # Write Hit, already has exclusive ownership. No message.
            self.memory_value = value
            self._print_status(f"{p_id} has a Write Hit (in M). No coherence message is sent.")
        elif state == 'E':
            # Write Hit, has exclusive ownership. No message.
            self.cache_states[p_id] = 'M'
            self.memory_value = value
            self._print_status(f"{p_id} has a Write Hit (in E). No message. State changes to M.")
        elif state == 'S':
            # Write Hit, but needs to invalidate others. Message sent.
            self.message_count += 1
            self.message_log.append(1)
            self.cache_states[p_id] = 'M'
            # Invalidate all other shared copies
            for p in self.processors:
                if p != p_id and self.cache_states[p] == 'S':
                    self.cache_states[p] = 'I'
            self.memory_value = value
            desc = (f"{p_id} has a Write Hit (in S). Sends 'Invalidate' message (1) "
                    f"to acquire exclusive ownership. State becomes M.")
            self._print_status(desc)
        elif state == 'I':
            # Write Miss. Message sent.
            self.message_count += 1
            self.message_log.append(1)
            # Invalidate all other copies
            for p in self.processors:
                if p != p_id:
                    self.cache_states[p] = 'I'
            self.cache_states[p_id] = 'M'
            self.memory_value = value
            desc = (f"{p_id} has a Write Miss. Sends 'Read For Ownership' message (1) "
                    f"to get data and invalidate others. State becomes M.")
            self._print_status(desc)
        print("-" * 25)
        
    def get_final_result(self):
        """Prints the final summary of the simulation."""
        print("--- Final Result ---")
        equation = " + ".join(map(str, self.message_log))
        print(f"Total messages = {equation} = {self.message_count}")


def solve():
    """
    Main function to run the MESI simulation based on the problem description.
    """
    processors = ['P1', 'P2', 'P3', 'P4']
    sim = MESISimulator(processors)

    # 1. P1 reads X
    sim.read('P1')
    
    # 2. P2 reads X
    sim.read('P2')

    # 3. P1 writes X = 1
    sim.write('P1', 1)

    # 4. P3 reads X
    sim.read('P3')

    # 5. P2 writes X = 2
    sim.write('P2', 2)

    # 6. P4 reads X
    sim.read('P4')

    # 7. P1 reads X
    sim.read('P1')

    # Print the final calculation
    sim.get_final_result()

solve()
<<<7>>>