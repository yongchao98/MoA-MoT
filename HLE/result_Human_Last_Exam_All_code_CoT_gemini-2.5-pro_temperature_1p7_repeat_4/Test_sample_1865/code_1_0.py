def solve_mesi_messages():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # Processor cache states for variable X, initialized to Invalid
    states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    message_count = 0
    message_log = []

    def print_status(operation, message_type, num_messages=1):
        nonlocal message_count
        message_count += num_messages
        message_log.append(str(num_messages))
        print(f"Operation: {operation}")
        print(f"  - Bus Message: {message_type}")
        print(f"  - New States: P1={states['P1']}, P2={states['P2']}, P3={states['P3']}, P4={states['P4']}")
        print(f"  - Messages for this step: {num_messages}")
        print(f"  - Cumulative Messages: {message_count}\n")

    print("--- Initial State ---")
    print(f"States: P1={states['P1']}, P2={states['P2']}, P3={states['P3']}, P4={states['P4']}")
    print(f"Cumulative Messages: {message_count}\n")

    # 1. P1 reads X
    # Read miss. P1 issues BusRd. No other cache has it. P1 -> E.
    states['P1'] = 'E'
    print_status("1. P1 reads X", "BusRd")

    # 2. P2 reads X
    # Read miss. P2 issues BusRd. P1 snoops, has E state. P1 -> S, P2 -> S.
    states['P1'] = 'S'
    states['P2'] = 'S'
    print_status("2. P2 reads X", "BusRd")

    # 3. P1 writes X = 1
    # Write hit on S state. P1 issues BusUpgr (Invalidate). P2 -> I. P1 -> M.
    states['P1'] = 'M'
    states['P2'] = 'I'
    print_status("3. P1 writes X = 1", "BusUpgr")

    # 4. P3 reads X
    # Read miss. P3 issues BusRd. P1 snoops (M state), writes back, supplies data.
    # P1 -> S, P3 -> S.
    states['P1'] = 'S'
    states['P3'] = 'S'
    print_status("4. P3 reads X", "BusRd")

    # 5. P2 writes X = 2
    # Write miss. P2 issues BusRdX. Sharers (P1, P3) are invalidated. P2 -> M.
    states['P1'] = 'I'
    states['P2'] = 'M'
    states['P3'] = 'I'
    print_status("5. P2 writes X = 2", "BusRdX")

    # 6. P4 reads X
    # Read miss. P4 issues BusRd. P2 snoops (M state), writes back, supplies data.
    # P2 -> S, P4 -> S.
    states['P2'] = 'S'
    states['P4'] = 'S'
    print_status("6. P4 reads X", "BusRd")

    # 7. P1 reads X
    # Read miss. P1 issues BusRd. Sharers (P2, P4) exist. One supplies data. P1 -> S.
    states['P1'] = 'S'
    print_status("7. P1 reads X", "BusRd")
    
    # Final result
    equation = " + ".join(message_log)
    print("--- Final Calculation ---")
    print(f"Total messages = {equation} = {message_count}")
    
    print(f"\n<<<{message_count}>>>")


solve_mesi_messages()