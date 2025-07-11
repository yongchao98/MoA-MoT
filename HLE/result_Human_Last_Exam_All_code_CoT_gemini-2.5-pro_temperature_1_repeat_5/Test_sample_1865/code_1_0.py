def solve_mesi_coherence():
    """
    Calculates the number of cache coherence messages for a given sequence of operations
    in a 4-processor system using the MESI protocol.
    """
    # Cache states for variable X in each processor {Processor: State}
    # Initial state: I (Invalid) for all
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    message_count = 0
    message_components = []

    def print_status(step, operation, message_type, messages_added):
        nonlocal message_count
        message_count += messages_added
        if messages_added > 0:
            message_components.append(str(messages_added))
            
        print(f"--- Step {step}: {operation} ---")
        if message_type:
            print(f"Message: {message_type} ({messages_added} message{'s' if messages_added > 1 else ''})")
        else:
            print("Message: None (Local cache operation)")
        print(f"Cache States: P1={caches['P1']}, P2={caches['P2']}, P3={caches['P3']}, P4={caches['P4']}")
        print(f"Running Message Count: {message_count}\n")

    # --- Initial State ---
    print("--- Initial State ---")
    print(f"Cache States: P1={caches['P1']}, P2={caches['P2']}, P3={caches['P3']}, P4={caches['P4']}")
    print(f"Running Message Count: {message_count}\n")

    # 1. P1 reads X
    # P1 has a read miss. It sends a Read Miss (BusRd) request on the bus.
    # Memory responds. Since no other cache has X, P1's cache state becomes Exclusive (E).
    caches['P1'] = 'E'
    print_status(1, "P1 reads X", "Read Miss (BusRd)", 1)

    # 2. P2 reads X
    # P2 has a read miss. It sends a Read Miss (BusRd) request.
    # P1 snoops the bus, sees the request, and supplies the data.
    # Both P1 and P2's cache states become Shared (S).
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print_status(2, "P2 reads X", "Read Miss (BusRd)", 1)

    # 3. P1 writes X = 1
    # P1's cache is in Shared (S) state. To write, it must invalidate other copies.
    # P1 sends a Bus Upgrade (Invalidate) message on the bus.
    # P2's state becomes Invalid (I). P1's state becomes Modified (M).
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print_status(3, "P1 writes X = 1", "Bus Upgrade (Invalidate)", 1)

    # 4. P3 reads X
    # P3 has a read miss. It sends a Read Miss (BusRd) request.
    # P1 snoops, has the data in Modified (M) state. It writes X back to memory and sends data to P3.
    # P1's state becomes Shared (S). P3's state becomes Shared (S).
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print_status(4, "P3 reads X", "Read Miss (BusRd)", 1)

    # 5. P2 writes X = 2
    # P2's cache is Invalid (I), so this is a write miss.
    # P2 sends a Read With Intent to Modify (BusRdX) message.
    # This message fetches the data and invalidates all other copies (P1 and P3).
    # P1 and P3 states become Invalid (I). P2's state becomes Modified (M).
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print_status(5, "P2 writes X = 2", "Read With Intent to Modify (BusRdX)", 1)

    # 6. P4 reads X
    # P4 has a read miss. It sends a Read Miss (BusRd) request.
    # P2 snoops, has the data in Modified (M) state. It writes X back to memory and sends data to P4.
    # P2's state becomes Shared (S). P4's state becomes Shared (S).
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print_status(6, "P4 reads X", "Read Miss (BusRd)", 1)
    
    # 7. P1 reads X
    # P1 has a read miss. It sends a Read Miss (BusRd) request.
    # P2 and P4 are in Shared (S) state. One of them (or memory) supplies the data.
    # P1's state becomes Shared (S). P2 and P4 remain Shared (S).
    caches['P1'] = 'S'
    print_status(7, "P1 reads X", "Read Miss (BusRd)", 1)
    
    # --- Final Summary ---
    print("--- Final Calculation ---")
    final_equation = " + ".join(message_components)
    print(f"Each of the 7 operations required one bus transaction, which counts as one coherence message.")
    print(f"Total Messages = {final_equation} = {message_count}")

solve_mesi_coherence()
<<<7>>>