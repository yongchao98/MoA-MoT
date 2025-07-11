def mesi_simulation():
    """
    Simulates the MESI protocol for a sequence of operations and counts cache coherence messages.
    """
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    processors = ['P1', 'P2', 'P3', 'P4']
    message_log = []
    step_messages = []

    def print_state(step, operation_desc):
        print(f"Step {step}: {operation_desc}")
        print(f"  - Initial states: {caches}")

    def log_message(msg_type):
        message_log.append(msg_type)
        step_messages.append(1)

    print("Initial State:")
    print(f"  - Variable X value in memory = 0")
    print(f"  - All processor caches are 'Invalid' (I) for X.\n")
    print("--- Simulating Operations ---")

    # Step 1: P1 reads X
    step = 1
    op_desc = "P1 reads X"
    print_state(step, op_desc)
    # P1 has a read miss (I state)
    # P1 sends BusRd. Since no other cache has the data, P1's state becomes Exclusive.
    log_message("BusRd")
    caches['P1'] = 'E'
    print(f"  - P1 cache is 'Invalid', causing a read miss. P1 sends a 'BusRd' message.")
    print(f"  - No other cache has the data. P1 loads from memory and its state becomes 'Exclusive'.")
    print(f"  - Final states:   {caches}")
    print(f"  - Messages so far: {len(message_log)}\n")

    # Step 2: P2 reads X
    step = 2
    op_desc = "P2 reads X"
    print_state(step, op_desc)
    # P2 has a read miss (I state)
    # P2 sends BusRd. P1 snoops, provides data, and changes its state from E to S. P2 becomes S.
    log_message("BusRd")
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print(f"  - P2 cache is 'Invalid', causing a read miss. P2 sends a 'BusRd' message.")
    print(f"  - P1 snoops the bus, provides the data from its 'Exclusive' line.")
    print(f"  - P1 state changes E -> S. P2 state becomes 'Shared'.")
    print(f"  - Final states:   {caches}")
    print(f"  - Messages so far: {len(message_log)}\n")

    # Step 3: P1 writes X = 1
    step = 3
    op_desc = "P1 writes X = 1"
    print_state(step, op_desc)
    # P1 is in S state. It must send a message to invalidate other sharers.
    log_message("BusUpgr")
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print(f"  - P1 is 'Shared', so a write requires invalidating other copies. P1 sends 'BusUpgr'.")
    print(f"  - P2 receives the invalidation and changes its state S -> I.")
    print(f"  - P1's state becomes 'Modified'.")
    print(f"  - Final states:   {caches}")
    print(f"  - Messages so far: {len(message_log)}\n")

    # Step 4: P3 reads X
    step = 4
    op_desc = "P3 reads X"
    print_state(step, op_desc)
    # P3 has a read miss (I state)
    # P3 sends BusRd. P1 snoops, provides the modified data, and changes state M -> S. P3 becomes S.
    log_message("BusRd")
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print(f"  - P3 cache is 'Invalid', causing a read miss. P3 sends 'BusRd'.")
    print(f"  - P1 snoops, sees its state is 'Modified', provides the data, and writes back to memory.")
    print(f"  - P1 state changes M -> S. P3 state becomes 'Shared'.")
    print(f"  - Final states:   {caches}")
    print(f"  - Messages so far: {len(message_log)}\n")

    # Step 5: P2 writes X = 2
    step = 5
    op_desc = "P2 writes X = 2"
    print_state(step, op_desc)
    # P2 has a write miss (I state)
    # P2 sends BusRdX. P1 and P3 invalidate (S -> I). P2 becomes M.
    log_message("BusRdX")
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print(f"  - P2 cache is 'Invalid', causing a write miss. P2 sends 'BusRdX'.")
    print(f"  - P1 and P3 receive the invalidation and change their state S -> I.")
    print(f"  - P2 loads the data and its state becomes 'Modified'.")
    print(f"  - Final states:   {caches}")
    print(f"  - Messages so far: {len(message_log)}\n")

    # Step 6: P4 reads X
    step = 6
    op_desc = "P4 reads X"
    print_state(step, op_desc)
    # P4 has a read miss (I state)
    # P4 sends BusRd. P2 snoops, provides data from its M line, and changes state M -> S. P4 becomes S.
    log_message("BusRd")
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print(f"  - P4 cache is 'Invalid', causing a read miss. P4 sends 'BusRd'.")
    print(f"  - P2 snoops ('Modified'), provides data, and writes back to memory.")
    print(f"  - P2 state changes M -> S. P4 state becomes 'Shared'.")
    print(f"  - Final states:   {caches}")
    print(f"  - Messages so far: {len(message_log)}\n")
    
    # Step 7: P1 reads X
    step = 7
    op_desc = "P1 reads X"
    print_state(step, op_desc)
    # P1 has a read miss (I state)
    # P1 sends BusRd. P2 and P4 are S. One provides the data. P1 becomes S.
    log_message("BusRd")
    caches['P1'] = 'S'
    print(f"  - P1 cache is 'Invalid', causing a read miss. P1 sends 'BusRd'.")
    print(f"  - P2/P4 are 'Shared' and one of them provides the data. Their states remain 'Shared'.")
    print(f"  - P1 state becomes 'Shared'.")
    print(f"  - Final states:   {caches}")
    print(f"  - Messages so far: {len(message_log)}\n")

    # Final Calculation
    print("--- Final Calculation ---")
    equation = " + ".join(map(str, step_messages))
    total_messages = len(message_log)
    print(f"The number of messages exchanged in each step are: {equation}.")
    print(f"Total cache coherence messages = {total_messages}")

if __name__ == '__main__':
    mesi_simulation()
    print("<<<7>>>")