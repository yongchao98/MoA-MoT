def solve_mesi_messages():
    """
    Traces the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # Initial states: I = Invalid
    cache_states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    message_count = 0
    message_log = []

    def log_message(reason):
        nonlocal message_count
        message_count += 1
        message_log.append(1)
        print(f"  - {reason} (Message +1)")

    def print_state(operation):
        print(f"\nOperation: {operation}")
        print(f"  - Current States: P1={cache_states['P1']}, P2={cache_states['P2']}, P3={cache_states['P3']}, P4={cache_states['P4']}")

    print("--- MESI Coherence Protocol Simulation ---")
    print(f"Initial States: P1={cache_states['P1']}, P2={cache_states['P2']}, P3={cache_states['P3']}, P4={cache_states['P4']}")
    print(f"Initial Message Count: {message_count}")

    # 1. P1 reads X
    print_state("1. P1 reads X")
    log_message("P1 has a Read Miss, issues a Read Miss on the bus.")
    cache_states['P1'] = 'E'
    print("  - State Change: P1 becomes Exclusive (E).")

    # 2. P2 reads X
    print_state("2. P2 reads X")
    log_message("P2 has a Read Miss, issues a Read Miss on the bus.")
    print("  - P1 snoops, provides data, and transitions to Shared.")
    cache_states['P1'] = 'S'
    cache_states['P2'] = 'S'
    print("  - State Change: P1 -> Shared (S), P2 -> Shared (S).")

    # 3. P1 writes X = 1
    print_state("3. P1 writes X = 1")
    log_message("P1 is Shared, issues an Invalidate message to upgrade.")
    cache_states['P2'] = 'I'
    cache_states['P1'] = 'M'
    print("  - State Change: P2 -> Invalid (I), P1 -> Modified (M).")

    # 4. P3 reads X
    print_state("4. P3 reads X")
    log_message("P3 has a Read Miss, issues a Read Miss on the bus.")
    print("  - P1 (Modified) snoops, provides data, and writes back to memory.")
    cache_states['P1'] = 'S'
    cache_states['P3'] = 'S'
    print("  - State Change: P1 -> Shared (S), P3 -> Shared (S).")

    # 5. P2 writes X = 2
    print_state("5. P2 writes X = 2")
    log_message("P2 has a Write Miss, issues a Read for Ownership (BusRdX).")
    print("  - P1 and P3 snoop and invalidate their copies.")
    cache_states['P1'] = 'I'
    cache_states['P3'] = 'I'
    cache_states['P2'] = 'M'
    print("  - State Change: P1 -> Invalid (I), P3 -> Invalid (I), P2 -> Modified (M).")

    # 6. P4 reads X
    print_state("6. P4 reads X")
    log_message("P4 has a Read Miss, issues a Read Miss on the bus.")
    print("  - P2 (Modified) snoops, provides data, and writes back to memory.")
    cache_states['P2'] = 'S'
    cache_states['P4'] = 'S'
    print("  - State Change: P2 -> Shared (S), P4 -> Shared (S).")

    # 7. P1 reads X
    print_state("7. P1 reads X")
    log_message("P1 has a Read Miss, issues a Read Miss on the bus.")
    print("  - P2 or P4 (Shared) snoops and provides the data.")
    cache_states['P1'] = 'S'
    print("  - State Change: P1 -> Shared (S).")

    print("\n--- Final Calculation ---")
    equation = " + ".join(map(str, message_log))
    print(f"The total number of messages is the sum of messages from each step.")
    print(f"Final equation: {equation} = {message_count}")
    print(f"Total cache coherence messages: {message_count}")

solve_mesi_messages()
<<<7>>>