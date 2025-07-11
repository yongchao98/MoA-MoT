import sys
from io import StringIO

def solve_mesi_trace():
    """
    Simulates the MESI protocol for a given sequence of operations
    and calculates the number of coherence messages.
    """

    # Initial state
    processors = ['P1', 'P2', 'P3', 'P4']
    caches = {p: {'state': 'I', 'value': None} for p in processors}
    memory_value = 0
    total_messages = 0
    messages_per_step = []

    # Helper function to print the current state
    def print_state(operation, message_generated):
        nonlocal total_messages
        nonlocal messages_per_step
        
        current_messages = 0
        if message_generated:
            total_messages += 1
            current_messages = 1
            messages_per_step.append(current_messages)

        states = ", ".join([f"{p}: {caches[p]['state']}" for p in processors])
        print(f"Step: {operation:<15} | States: ({states}) | Message: {message_generated} ({current_messages} sent)")

    # Sequence of operations
    print("--- MESI Simulation Trace ---")
    print(f"Initial State          | States: (P1: I, P2: I, P3: I, P4: I) | Message: False (0 sent)")

    # 1. P1 reads X
    # P1 has a read miss, issues BusRd. No other cache has X.
    # P1 gets data from memory and transitions to Exclusive (E).
    caches['P1']['state'] = 'E'
    caches['P1']['value'] = memory_value
    print_state("P1 reads X", True)

    # 2. P2 reads X
    # P2 has a read miss, issues BusRd. P1 snoops and has X in state E.
    # P1 supplies data to P2 and transitions to Shared (S). P2 transitions to Shared (S).
    caches['P1']['state'] = 'S'
    caches['P2']['state'] = 'S'
    caches['P2']['value'] = caches['P1']['value']
    print_state("P2 reads X", True)

    # 3. P1 writes X = 1
    # P1 has X in state S. To write, it must upgrade. Issues Invalidate message.
    # P2 snoops, sees the Invalidate, and sets its copy to Invalid (I).
    # P1 transitions to Modified (M).
    caches['P2']['state'] = 'I'
    caches['P1']['state'] = 'M'
    caches['P1']['value'] = 1
    print_state("P1 writes X=1", True)

    # 4. P3 reads X
    # P3 has a read miss, issues BusRd. P1 snoops and has X in state M.
    # P1 supplies the data (value 1) to P3 and main memory.
    # P1 transitions to Shared (S). P3 transitions to Shared (S).
    memory_value = caches['P1']['value'] # P1 flushes to memory
    caches['P1']['state'] = 'S'
    caches['P3']['state'] = 'S'
    caches['P3']['value'] = caches['P1']['value']
    print_state("P3 reads X", True)

    # 5. P2 writes X = 2
    # P2 has a write miss (state I), issues BusRdX.
    # P1 and P3 snoop, see the BusRdX, and invalidate their copies.
    # P2 fetches the block, writes the new value, and transitions to Modified (M).
    caches['P1']['state'] = 'I'
    caches['P3']['state'] = 'I'
    caches['P2']['state'] = 'M'
    caches['P2']['value'] = 2
    print_state("P2 writes X=2", True)

    # 6. P4 reads X
    # P4 has a read miss, issues BusRd. P2 snoops and has X in state M.
    # P2 supplies the data (value 2) to P4 and main memory.
    # P2 transitions to Shared (S). P4 transitions to Shared (S).
    memory_value = caches['P2']['value'] # P2 flushes to memory
    caches['P2']['state'] = 'S'
    caches['P4']['state'] = 'S'
    caches['P4']['value'] = caches['P2']['value']
    print_state("P4 reads X", True)

    # 7. P1 reads X
    # P1 has a read miss (state I), issues BusRd.
    # P2 and P4 have the data in state S. One of them supplies the data.
    # P1 transitions to Shared (S). P2 and P4 remain Shared (S).
    caches['P1']['state'] = 'S'
    caches['P1']['value'] = caches['P2']['value']
    print_state("P1 reads X", True)

    print("\n--- Final Calculation ---")
    equation = " + ".join(map(str, messages_per_step))
    print(f"The total number of cache coherence messages is the sum of messages at each step:")
    print(f"{equation} = {total_messages}")

    # The final required output format
    return total_messages

if __name__ == '__main__':
    result = solve_mesi_trace()
    print(f'<<<{result}>>>')
