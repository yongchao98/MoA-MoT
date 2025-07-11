import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

class MesiSimulator:
    """
    Simulates the MESI protocol for a given sequence of operations.
    Counts the number of cache coherence messages.
    """
    def __init__(self, num_processors):
        # Initial state: All caches are Invalid (I) for the variable X
        self.processors = {f"P{i+1}": "I" for i in range(num_processors)}
        self.message_count = 0
        self.message_log = []
        self.step_counter = 0

    def log_and_update(self, operation_desc, bus_op, messages, message_count):
        """Logs the details of a single step and updates the total message count."""
        self.step_counter += 1
        print(f"Step {self.step_counter}: {operation_desc}")
        print(f"  Initial States: {self.initial_states}")
        print(f"  Bus Operation: {bus_op}")
        print(f"  Coherence Messages: {messages} ({message_count} message{'s' if message_count != 1 else ''})")
        self.message_count += message_count
        self.message_log.append(message_count)
        print(f"  Final States:   {self.processors}")
        print(f"  Cumulative Messages: {self.message_count}\n")

    def read(self, p_id):
        """Simulates a processor read operation."""
        self.initial_states = self.processors.copy()
        operation_desc = f"P{p_id[-1]} reads X"

        # Cache hit: No messages if state is M, E, or S
        if self.processors[p_id] in ["M", "E", "S"]:
            self.log_and_update(operation_desc, "Cache Hit", "None", 0)
            return

        # Cache miss: A bus operation is required
        other_states = {p: s for p, s in self.processors.items() if p != p_id}
        
        # Case 1: Another cache has the data in Modified (M) state
        if "M" in other_states.values():
            holder = [p for p, s in other_states.items() if s == "M"][0]
            self.processors[holder] = "S"
            self.processors[p_id] = "S"
            self.log_and_update(operation_desc, "Read Miss", "Request from P" + p_id[-1] + ", Data Response from " + holder, 2)
        # Case 2: Another cache has the data in Exclusive (E) state
        elif "E" in other_states.values():
            holder = [p for p, s in other_states.items() if s == "E"][0]
            self.processors[holder] = "S"
            self.processors[p_id] = "S"
            self.log_and_update(operation_desc, "Read Miss", "Request from P" + p_id[-1] + ", Data Response from " + holder, 2)
        # Case 3: Other caches have the data in Shared (S) state
        elif "S" in other_states.values():
            holder = [p for p, s in other_states.items() if s == "S"][0]
            self.processors[p_id] = "S"
            self.log_and_update(operation_desc, "Read Miss", "Request from P" + p_id[-1] + ", Data Response from " + holder, 2)
        # Case 4: No other cache has the data (fetch from memory)
        else:
            self.processors[p_id] = "E"
            self.log_and_update(operation_desc, "Read Miss", "Request from P" + p_id[-1] + ", Data Response from Memory", 2)

    def write(self, p_id, value):
        """Simulates a processor write operation."""
        self.initial_states = self.processors.copy()
        operation_desc = f"P{p_id[-1]} writes X = {value}"

        current_state = self.processors[p_id]

        # Case 1: Write hit (state is M or E). No bus messages.
        if current_state in ["M", "E"]:
            self.processors[p_id] = "M"
            self.log_and_update(operation_desc, "Cache Hit (No broadcast needed)", "None", 0)
            return
            
        # Case 2: Write hit (state is S). Needs to invalidate other copies.
        if current_state == "S":
            for p in self.processors:
                if p != p_id and self.processors[p] == "S":
                    self.processors[p] = "I"
            self.processors[p_id] = "M"
            self.log_and_update(operation_desc, "Upgrade", "Invalidate broadcast from P" + p_id[-1], 1)
            
        # Case 3: Write miss (state is I). Needs to get data and invalidate others.
        elif current_state == "I":
            # Invalidate all other caches
            for p in self.processors:
                if p != p_id:
                    self.processors[p] = "I"
            self.processors[p_id] = "M"
            self.log_and_update(operation_desc, "Read For Ownership (RFO)", "RFO Request from P" + p_id[-1] + ", Data Response", 2)

    def run_simulation(self):
        """Executes the full sequence of operations."""
        print("--- MESI Cache Coherence Simulation ---\n")
        print(f"Initial States: {self.processors}\n")

        # 1. P1 reads X
        self.read("P1")
        # 2. P2 reads X
        self.read("P2")
        # 3. P1 writes X = 1
        self.write("P1", 1)
        # 4. P3 reads X
        self.read("P3")
        # 5. P2 writes X = 2
        self.write("P2", 2)
        # 6. P4 reads X
        self.read("P4")
        # 7. P1 reads X
        self.read("P1")

        print("--- Simulation Complete ---\n")
        print(f"Final Cache States: {self.processors}")
        
        # Print the final calculation
        equation = " + ".join(map(str, self.message_log))
        print(f"Total Message Calculation: {equation} = {self.message_count}")


# Main execution
if __name__ == "__main__":
    simulator = MesiSimulator(num_processors=4)
    simulator.run_simulation()
    
    # Restore stdout and print the captured output
    sys.stdout = old_stdout
    output = captured_output.getvalue()
    print(output)
