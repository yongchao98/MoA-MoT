def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # MESI states: Modified, Exclusive, Shared, Invalid
    processors = ['P1', 'P2', 'P3', 'P4']
    caches = {p: 'I' for p in processors}
    message_count = 0
    messages_per_step = []

    def print_state(operation_description):
        print(f"\nOperation: {operation_description}")
        print(f"  Cache States: {caches}")
        print(f"  Messages so far: {message_count}")

    print("Initial State: X is not in any cache, value in memory is 0.")
    print(f"  Cache States: {caches}")
    print(f"  Message Count: {message_count}")

    # 1. P1 reads X
    # P1 has a read miss. It sends a Read Miss (BusRd) message.
    # Since no other cache has X, memory responds. P1's cache state becomes Exclusive (E).
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'E'
    print_state("1. P1 reads X")
    print("  -> P1 has a read miss, sends a 'Read Miss' message.")
    print("  -> P1's state for X becomes Exclusive (E).")

    # 2. P2 reads X
    # P2 has a read miss. It sends a Read Miss (BusRd) message.
    # P1 snoops the bus, sees the request, and provides the data.
    # Both P1 and P2's cache states become Shared (S).
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print_state("2. P2 reads X")
    print("  -> P2 has a read miss, sends a 'Read Miss' message.")
    print("  -> P1 provides the data. P1 and P2 states for X become Shared (S).")

    # 3. P1 writes X = 1
    # P1's cache has X in Shared (S) state. To write, it must have exclusive ownership.
    # P1 sends an Invalidate (BusUpgr) message.
    # P2 receives the invalidate message and sets its state for X to Invalid (I).
    # P1's state becomes Modified (M).
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print_state("3. P1 writes X = 1")
    print("  -> P1 needs to write to a Shared line, sends an 'Invalidate' message.")
    print("  -> P2's copy is invalidated. P1's state for X becomes Modified (M).")

    # 4. P3 reads X
    # P3 has a read miss. It sends a Read Miss (BusRd) message.
    # P1 has the data in Modified (M) state. It provides the data to P3 and writes it back to memory.
    # Both P1 and P3's cache states become Shared (S).
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print_state("4. P3 reads X")
    print("  -> P3 has a read miss, sends a 'Read Miss' message.")
    print("  -> P1 provides the data. P1 and P3 states for X become Shared (S).")

    # 5. P2 writes X = 2
    # P2 has a write miss (its state is I). It sends a Write Miss (BusRdX) message.
    # This message invalidates all other shared copies.
    # P1 and P3 set their states to Invalid (I).
    # P2's state becomes Modified (M).
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print_state("5. P2 writes X = 2")
    print("  -> P2 has a write miss, sends a 'Write Miss' (BusRdX) message.")
    print("  -> P1 and P3 copies are invalidated. P2's state for X becomes Modified (M).")

    # 6. P4 reads X
    # P4 has a read miss. It sends a Read Miss (BusRd) message.
    # P2 has the data in Modified (M) state. It provides the data to P4 and writes it back to memory.
    # Both P2 and P4's cache states become Shared (S).
    message_count += 1
    messages_per_step.append(1)
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print_state("6. P4 reads X")
    print("  -> P4 has a read miss, sends a 'Read Miss' message.")
    print("  -> P2 provides the data. P2 and P4 states for X become Shared (S).")

    # 7. P1 reads X
    # P1 has a read miss (its state is I). It sends a Read Miss (BusRd) message.
    # P2 and P4 have the data in Shared (S) state. One of them (or memory) provides the data.
    # P1's state becomes Shared (S).
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'S'
    print_state("7. P1 reads X")
    print("  -> P1 has a read miss, sends a 'Read Miss' message.")
    print("  -> Data is provided by a shared cache. P1's state for X becomes Shared (S).")

    # Final Calculation
    print("\n" + "="*30)
    print("Final Calculation Summary")
    print("="*30)
    print("Each operation that required bus access generated one coherence message:")
    equation = " + ".join(map(str, messages_per_step))
    print(f"Total Messages = {equation} = {message_count}")

solve_mesi_coherence()
print("\n<<<7>>>")