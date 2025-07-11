def solve_mesi_messages():
    """
    Simulates the MESI protocol for a sequence of operations and counts the
    cache coherence messages.
    """
    # Cache states: I (Invalid), S (Shared), E (Exclusive), M (Modified)
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    message_count = 0
    message_log = []

    def print_status(op_num, description, message_type):
        nonlocal message_count
        message_count += 1
        message_log.append(1)
        print(f"Step {op_num}: {description}")
        print(f"  - Action: This operation causes a {message_type}.")
        print(f"  - Message: 1 bus message is generated.")
        print(f"  - New States: {caches}")
        print("-" * 20)

    print("Initial State: X is not in any cache. All caches are 'Invalid'. Memory(X)=0\n")

    # 1. P1 reads X
    # P1 has a read miss. It issues a Read Miss on the bus.
    # No other cache has the data, so memory provides it.
    # P1's cache state for X becomes Exclusive (E).
    caches['P1'] = 'E'
    print_status(1, "P1 reads X", "Read Miss")

    # 2. P2 reads X
    # P2 has a read miss. It issues a Read Miss on the bus.
    # P1 snoops the bus and has the data in state E.
    # P1 provides the data to P2 (cache-to-cache transfer) and changes its state to Shared (S).
    # P2's cache state becomes Shared (S).
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print_status(2, "P2 reads X", "Read Miss")

    # 3. P1 writes X = 1
    # P1 wants to write, but its state is S. It must gain exclusive ownership.
    # P1 issues an Invalidate message on the bus.
    # P2 snoops the bus and invalidates its copy of X (S -> I).
    # P1's state changes to Modified (M).
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print_status(3, "P1 writes X = 1", "Write Miss (Invalidate)")

    # 4. P3 reads X
    # P3 has a read miss. It issues a Read Miss on the bus.
    # P1 snoops and has the data in state M.
    # P1 provides the data to P3, writes the data back to main memory, and changes its state to S.
    # P3's state becomes S.
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print_status(4, "P3 reads X", "Read Miss")

    # 5. P2 writes X = 2
    # P2 has a write miss (its state is I). It needs to get the data and ownership.
    # P2 issues a Read For Ownership (RFO) message (which is also an invalidate).
    # P1 and P3 snoop and invalidate their copies (S -> I).
    # One cache (or memory) provides the data to P2. P2's state becomes M.
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print_status(5, "P2 writes X = 2", "Write Miss (Read For Ownership)")

    # 6. P4 reads X
    # P4 has a read miss. It issues a Read Miss on the bus.
    # P2 snoops and has the data in state M.
    # P2 provides the data to P4, writes back to memory, and changes its state to S.
    # P4's state becomes S.
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print_status(6, "P4 reads X", "Read Miss")

    # 7. P1 reads X
    # P1 has a read miss (its state is I). It issues a Read Miss on the bus.
    # P2 and P4 have the data in state S.
    # One of them (e.g., P2) provides the data to P1 via cache-to-cache transfer.
    # P1's state becomes S. The states of P2 and P4 remain S.
    caches['P1'] = 'S'
    print_status(7, "P1 reads X", "Read Miss")

    # Final Count
    equation = " + ".join(map(str, message_log))
    print(f"\nTotal messages = (Msg from step 1) + (Msg from step 2) + ... + (Msg from step 7)")
    print(f"Final Calculation: {equation} = {message_count}")
    print(f"\nTotal number of cache coherence messages exchanged is {message_count}.")

solve_mesi_messages()
<<<7>>>