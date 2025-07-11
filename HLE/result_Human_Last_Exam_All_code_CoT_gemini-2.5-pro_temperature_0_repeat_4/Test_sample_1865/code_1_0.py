def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a given sequence of operations
    and calculates the number of cache coherence messages.
    """
    # --- Initialization ---
    # Cache states for variable X in each processor
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    # Value of X in main memory
    memory_value = 0
    # Value of X in each processor's cache
    cache_values = {'P1': None, 'P2': None, 'P3': None, 'P4': None}
    
    message_count = 0
    messages_per_step = []

    def print_state(step_description, details, step_messages):
        """Helper function to print the state at each step."""
        print(f"--- {step_description} ---")
        print(details)
        print(f"Messages for this step: {step_messages}")
        print(f"Cache States: P1={caches['P1']}, P2={caches['P2']}, P3={caches['P3']}, P4={caches['P4']}")
        print(f"Total Messages so far: {message_count}\n")

    print("### MESI Cache Coherence Simulation ###\n")
    print("Initial State: X is not in any cache. Memory value of X = 0.\n")

    # --- Sequence of Operations ---

    # 1. P1 reads X
    step_messages = 1
    message_count += step_messages
    messages_per_step.append(step_messages)
    # P1 has a read miss. It issues a Bus Read (BusRd).
    # Memory responds. P1's state becomes Exclusive (E) as it's the only cache with the data.
    caches['P1'] = 'E'
    cache_values['P1'] = memory_value
    print_state("Step 1: P1 reads X",
                "P1 has a cache miss and issues a Bus Read (BusRd).\nMemory provides the data. P1's cache state for X becomes Exclusive (E).",
                step_messages)

    # 2. P2 reads X
    step_messages = 1
    message_count += step_messages
    messages_per_step.append(step_messages)
    # P2 has a read miss. It issues a Bus Read (BusRd).
    # P1 snoops the bus, sees the request, and provides the data (cache-to-cache transfer).
    # P1's state changes from E to Shared (S). P2's state becomes S.
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    cache_values['P2'] = cache_values['P1']
    print_state("Step 2: P2 reads X",
                "P2 has a cache miss and issues a Bus Read (BusRd).\nP1 snoops, provides the data, and changes its state from E to Shared (S).\nP2's cache state for X becomes Shared (S).",
                step_messages)

    # 3. P1 writes X = 1
    step_messages = 1
    message_count += step_messages
    messages_per_step.append(step_messages)
    # P1 has X in state S. To write, it must invalidate other copies.
    # P1 issues an Invalidate message (or BusUpgrade).
    # P2 snoops and invalidates its copy (S -> I).
    # P1's state changes from S to Modified (M).
    caches['P2'] = 'I'
    caches['P1'] = 'M'
    cache_values['P1'] = 1
    print_state("Step 3: P1 writes X = 1",
                "P1 has X in state S and needs to write. It issues an Invalidate message.\nP2 invalidates its copy (S -> I). P1's state becomes Modified (M).",
                step_messages)

    # 4. P3 reads X
    step_messages = 1
    message_count += step_messages
    messages_per_step.append(step_messages)
    # P3 has a read miss. It issues a Bus Read (BusRd).
    # P1 snoops, sees it has X in state M. It provides the data and writes back to memory.
    # P1's state changes from M to S. P3's state becomes S.
    memory_value = cache_values['P1']
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    cache_values['P3'] = cache_values['P1']
    print_state("Step 4: P3 reads X",
                "P3 has a cache miss and issues a Bus Read (BusRd).\nP1 snoops, provides the data (and writes back to memory), and changes its state from M to S.\nP3's cache state for X becomes Shared (S).",
                step_messages)

    # 5. P2 writes X = 2
    step_messages = 1
    message_count += step_messages
    messages_per_step.append(step_messages)
    # P2 has a write miss (state I). It issues a Bus Read Exclusive (BusRdX).
    # This message serves to get the data and invalidate other copies.
    # P1 and P3 snoop and invalidate their copies (S -> I).
    # P2's state becomes M.
    caches['P1'] = 'I'
    caches['P3'] = 'I'
    caches['P2'] = 'M'
    cache_values['P2'] = 2
    print_state("Step 5: P2 writes X = 2",
                "P2 has a write miss and issues a Bus Read Exclusive (BusRdX).\nP1 and P3 invalidate their copies (S -> I).\nP2's cache state for X becomes Modified (M).",
                step_messages)

    # 6. P4 reads X
    step_messages = 1
    message_count += step_messages
    messages_per_step.append(step_messages)
    # P4 has a read miss. It issues a Bus Read (BusRd).
    # P2 snoops, sees it has X in state M. It provides the data and writes back to memory.
    # P2's state changes from M to S. P4's state becomes S.
    memory_value = cache_values['P2']
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    cache_values['P4'] = cache_values['P2']
    print_state("Step 6: P4 reads X",
                "P4 has a cache miss and issues a Bus Read (BusRd).\nP2 snoops, provides the data (and writes back to memory), and changes its state from M to S.\nP4's cache state for X becomes Shared (S).",
                step_messages)

    # 7. P1 reads X
    step_messages = 1
    message_count += step_messages
    messages_per_step.append(step_messages)
    # P1 has a read miss (state I). It issues a Bus Read (BusRd).
    # P2 and P4 have the data in state S. One of them provides it.
    # P1's state becomes S.
    caches['P1'] = 'S'
    cache_values['P1'] = cache_values['P2']
    print_state("Step 7: P1 reads X",
                "P1 has a cache miss and issues a Bus Read (BusRd).\nA current sharer (P2 or P4) provides the data.\nP1's cache state for X becomes Shared (S).",
                step_messages)

    # --- Final Calculation ---
    print("--- Final Calculation ---")
    equation = " + ".join(map(str, messages_per_step))
    print(f"The total number of cache coherence messages is the sum of messages from each step.")
    print(f"Total messages = {equation} = {message_count}")

if __name__ == '__main__':
    solve_mesi_coherence()
    print("\n<<<7>>>")