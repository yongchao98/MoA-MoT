def simulate_mesi():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    processors = ['P1', 'P2', 'P3', 'P4']
    message_count = 0
    message_log = []

    def print_state(op_num, description):
        states = ", ".join([f"{p}: {caches[p]}" for p in processors])
        print(f"Step {op_num}: {description}")
        print(f"  -> States: [{states}] | Messages: {message_count}\n")

    print("Initial State:")
    print_state(0, "All caches are Invalid.")

    # 1. P1 reads X
    message_count += 1
    message_log.append(1)
    # P1 has a read miss, issues BusRd. No other cache has X. P1 state -> E.
    caches['P1'] = 'E'
    print_state(1, "P1 reads X. Issues BusRd, moves to Exclusive.")

    # 2. P2 reads X
    message_count += 1
    message_log.append(1)
    # P2 has a read miss, issues BusRd. P1 snoops, provides data.
    # P1 state E -> S. P2 state -> S.
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print_state(2, "P2 reads X. Issues BusRd. P1 provides data. P1, P2 move to Shared.")

    # 3. P1 writes X = 1
    message_count += 1
    message_log.append(1)
    # P1 has write hit on Shared line. Issues Invalidate message.
    # P2 state S -> I. P1 state S -> M.
    caches['P2'] = 'I'
    caches['P1'] = 'M'
    print_state(3, "P1 writes X. Issues Invalidate. P1 moves to Modified, P2 invalidates.")

    # 4. P3 reads X
    message_count += 1
    message_log.append(1)
    # P3 has read miss, issues BusRd. P1 snoops, provides data.
    # P1 state M -> S. P3 state -> S.
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print_state(4, "P3 reads X. Issues BusRd. P1 provides data. P1, P3 move to Shared.")

    # 5. P2 writes X = 2
    message_count += 1
    message_log.append(1)
    # P2 has write miss. Issues BusRdX (Read for Ownership).
    # P1 and P3 snoop, invalidate their copies.
    # P1 state S -> I. P3 state S -> I. P2 state -> M.
    caches['P1'] = 'I'
    caches['P3'] = 'I'
    caches['P2'] = 'M'
    print_state(5, "P2 writes X. Issues BusRdX. P2 moves to Modified, P1 and P3 invalidate.")

    # 6. P4 reads X
    message_count += 1
    message_log.append(1)
    # P4 has read miss, issues BusRd. P2 snoops, provides data.
    # P2 state M -> S. P4 state -> S.
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print_state(6, "P4 reads X. Issues BusRd. P2 provides data. P2, P4 move to Shared.")

    # 7. P1 reads X
    message_count += 1
    message_log.append(1)
    # P1 has read miss, issues BusRd. P2 or P4 provides data.
    # P1 state -> S. P2 and P4 remain S.
    caches['P1'] = 'S'
    print_state(7, "P1 reads X. Issues BusRd. P2/P4 provides data. P1 moves to Shared.")
    
    # Final Count
    equation = " + ".join(map(str, message_log))
    print("--------------------------------------------------")
    print(f"Total number of cache coherence messages:")
    print(f"{equation} = {message_count}")
    print("--------------------------------------------------")

simulate_mesi()