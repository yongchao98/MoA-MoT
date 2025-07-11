def solve_mesi_simulation():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # Initialization
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    memory_value = 0
    message_count = 0
    message_log = []

    def print_state(operation_description):
        print(f"--- After {operation_description} ---")
        state_str = ", ".join([f"{p}: {s}" for p, s in caches.items()])
        print(f"Cache States: {state_str}")
        print(f"Messages so far: {message_count}\n")

    print("--- Initial State ---")
    print(f"Cache States: P1: I, P2: I, P3: I, P4: I")
    print(f"Memory Value of X: {memory_value}\n")

    # 1. P1 reads X
    # P1 has a Read Miss. It sends a 'Read' message on the bus.
    # No other cache has X, so P1 gets it from memory and enters the 'Exclusive' (E) state.
    message_count += 1
    message_log.append("1 (P1 Read Miss -> BusRd)")
    caches['P1'] = 'E'
    print_state("1. P1 reads X")

    # 2. P2 reads X
    # P2 has a Read Miss. It sends a 'Read' message on the bus.
    # P1 snoops and provides the data. P1's state changes from E -> 'Shared' (S).
    # P2 receives the data and enters the 'Shared' (S) state.
    message_count += 1
    message_log.append("1 (P2 Read Miss -> BusRd)")
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print_state("2. P2 reads X")

    # 3. P1 writes X = 1
    # P1 wants to write but is in state S. It must upgrade to 'Modified' (M).
    # P1 sends an 'Invalidate' message on the bus.
    # P2 snoops, sees the invalidate, and changes its state S -> 'Invalid' (I).
    # P1 changes its state S -> M and performs the write.
    message_count += 1
    message_log.append("1 (P1 Write Hit -> BusUpgr/Invalidate)")
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print_state("3. P1 writes X = 1")

    # 4. P3 reads X
    # P3 has a Read Miss. It sends a 'Read' message.
    # P1 snoops (state M), provides the data to P3, and writes the data back to memory.
    # P1's state changes M -> S. P3 enters state S.
    message_count += 1
    message_log.append("1 (P3 Read Miss -> BusRd)")
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print_state("4. P3 reads X")

    # 5. P2 writes X = 2
    # P2 has a Write Miss (state I).
    # P2 sends a 'Read for Ownership' (RFO) message (also called BusRdX).
    # P1 and P3 snoop, see the RFO, and invalidate their copies (S -> I).
    # One of them provides the data to P2, which P2 then overwrites. P2 enters state M.
    message_count += 1
    message_log.append("1 (P2 Write Miss -> BusRdX/RFO)")
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print_state("5. P2 writes X = 2")

    # 6. P4 reads X
    # P4 has a Read Miss. It sends a 'Read' message.
    # P2 snoops (state M), provides the data to P4, and writes the data back to memory.
    # P2's state changes M -> S. P4 enters state S.
    message_count += 1
    message_log.append("1 (P4 Read Miss -> BusRd)")
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print_state("6. P4 reads X")

    # 7. P1 reads X
    # P1 has a Read Miss (state I). It sends a 'Read' message.
    # P2 and P4 have the data in state S. One of them provides the data to P1.
    # P1 enters state S. P2 and P4 remain in S.
    message_count += 1
    message_log.append("1 (P1 Read Miss -> BusRd)")
    caches['P1'] = 'S'
    print_state("7. P1 reads X")

    # Final Calculation
    print("--- Final Calculation ---")
    equation = " + ".join(message_log)
    print(f"Total messages = {equation} = {message_count}")
    return message_count

if __name__ == '__main__':
    total_messages = solve_mesi_simulation()
    print(f"\n<<<7>>>")