def solve_mesi_coherence():
    """
    Analyzes a sequence of operations to calculate the number of MESI
    cache coherence messages.

    The function simulates the state changes of a shared variable 'X' in four
    processor caches (P1, P2, P3, P4) and counts the messages exchanged over
    the bus for each operation.

    Assumptions:
    - MESI protocol is used.
    - Cache-to-cache transfers are allowed and preferred over memory access.
    - A read/write request to the bus is one message.
    - A data response (from memory or another cache) is one message.
    - A write-back to memory is one message.
    - An invalidate signal (or a request like RFO that implies invalidation)
      is counted as one message that is broadcast on the bus.
    """
    
    # M: Modified, E: Exclusive, S: Shared, I: Invalid
    cache_states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    memory_value = 0
    total_messages = 0
    message_counts = []

    def print_state(step, operation, messages_this_step, explanation):
        nonlocal total_messages
        total_messages += messages_this_step
        message_counts.append(messages_this_step)
        
        print(f"--- Step {step}: {operation} ---")
        print(explanation)
        print(f"Messages this step: {messages_this_step}")
        print(f"Current Cache States: P1={cache_states['P1']}, P2={cache_states['P2']}, P3={cache_states['P3']}, P4={cache_states['P4']}")
        print(f"Running Total Messages: {total_messages}\n")

    # Initial State
    print("Initial State:")
    print(f"Cache States: P1={cache_states['P1']}, P2={cache_states['P2']}, P3={cache_states['P3']}, P4={cache_states['P4']}")
    print(f"Memory Value of X: {memory_value}")
    print(f"Total Messages: {total_messages}\n")

    # 1. P1 reads X
    # P1 has a cache miss.
    # - P1 sends a Read request to the bus (1 message).
    # - No other cache has the data. Memory responds with the data (1 message).
    # P1 is the only cache with the data, so its state becomes Exclusive.
    step = 1
    operation = "P1 reads X"
    messages = 2
    explanation = "P1 has a read miss. It sends a 'Read' request to the bus (1). Memory provides the data (1)."
    cache_states['P1'] = 'E'
    print_state(step, operation, messages, explanation)

    # 2. P2 reads X
    # P2 has a cache miss.
    # - P2 sends a Read request to the bus (1 message).
    # - P1 snoops the bus, sees the request, and has the data in state E.
    # - P1 provides the data to P2 (cache-to-cache transfer) (1 message).
    # P1's state changes from E to Shared. P2's state becomes Shared.
    step = 2
    operation = "P2 reads X"
    messages = 2
    explanation = "P2 has a read miss. It sends a 'Read' request to the bus (1). P1 snoops and supplies the data (1). Both caches are now in the 'Shared' state."
    cache_states['P1'] = 'S'
    cache_states['P2'] = 'S'
    print_state(step, operation, messages, explanation)
    
    # 3. P1 writes X = 1
    # P1 has the data in Shared state, which is read-only.
    # - P1 sends an 'Invalidate' or 'Upgrade' signal to the bus to gain exclusive ownership (1 message).
    # - P2 receives the signal and invalidates its copy (S -> I).
    # P1's state changes from S to Modified as it writes the new value.
    step = 3
    operation = "P1 writes X = 1"
    messages = 1
    explanation = "P1 must upgrade its state to 'Modified' to write. It sends an 'Invalidate' message on the bus (1). P2 invalidates its copy."
    cache_states['P1'] = 'M'
    cache_states['P2'] = 'I'
    memory_value = 1
    print_state(step, operation, messages, explanation)
    
    # 4. P3 reads X
    # P3 has a cache miss.
    # - P3 sends a Read request to the bus (1 message).
    # - P1 snoops, has the data in M state. It must intervene.
    # - P1 writes the modified data back to memory (1 message).
    # - P1 sends the data to P3 (1 message).
    # P1's state changes from M to S. P3's state becomes S.
    step = 4
    operation = "P3 reads X"
    messages = 3
    explanation = "P3 has a read miss. It sends 'Read' (1). P1 has the data in 'Modified', so it writes back to memory (1) and supplies the data to P3 (1)."
    cache_states['P1'] = 'S'
    cache_states['P3'] = 'S'
    print_state(step, operation, messages, explanation)

    # 5. P2 writes X = 2
    # P2's copy is Invalid (a write miss).
    # - P2 sends a 'Read For Ownership' (RFO) request to the bus (1 message). This request claims the line and invalidates other copies.
    # - P1 and P3 are Shared. One of them (e.g., P1) supplies the data to P2 (1 message).
    # - P1 and P3 see the RFO and invalidate their copies (S -> I).
    # P2's state becomes Modified.
    step = 5
    operation = "P2 writes X = 2"
    messages = 2
    explanation = "P2 has a write miss. It sends a 'Read For Ownership' (RFO) request (1). A shared cache (P1 or P3) supplies the data (1). Both P1 and P3 invalidate their copies."
    cache_states['P1'] = 'I'
    cache_states['P2'] = 'M'
    cache_states['P3'] = 'I'
    memory_value = 2
    print_state(step, operation, messages, explanation)

    # 6. P4 reads X
    # P4 has a cache miss.
    # - P4 sends a Read request to the bus (1 message).
    # - P2 snoops, has the data in M state.
    # - P2 writes the data back to memory (1 message).
    # - P2 sends the data to P4 (1 message).
    # P2's state changes from M to S. P4's state becomes S.
    step = 6
    operation = "P4 reads X"
    messages = 3
    explanation = "P4 has a read miss. It sends 'Read' (1). P2 has data in 'Modified', so it writes back to memory (1) and supplies data to P4 (1)."
    cache_states['P2'] = 'S'
    cache_states['P4'] = 'S'
    print_state(step, operation, messages, explanation)

    # 7. P1 reads X
    # P1's copy is Invalid.
    # - P1 sends a Read request to the bus (1 message).
    # - P2 and P4 have the data in Shared state. One of them (e.g., P2) supplies the data (1 message).
    # P1's state becomes S. Others remain S.
    step = 7
    operation = "P1 reads X"
    messages = 2
    explanation = "P1 has a read miss. It sends 'Read' (1). A shared cache (P2 or P4) supplies the data (1). P1's state becomes 'Shared'."
    cache_states['P1'] = 'S'
    print_state(step, operation, messages, explanation)

    # Final summary
    print("--- Final Calculation ---")
    equation = " + ".join(map(str, message_counts))
    print(f"Total Messages = {equation} = {total_messages}")
    print(f"<<<{total_messages}>>>")

if __name__ == '__main__':
    solve_mesi_coherence()