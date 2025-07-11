def solve_mesi_simulation():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # Initial State: [P1, P2, P3, P4]
    # 'I' for Invalid, 'S' for Shared, 'E' for Exclusive, 'M' for Modified
    cache_states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    messages = 0
    message_counts = []

    print("--- MESI Cache Coherence Simulation ---")
    print("Initial State:")
    print(f"  Cache States: {cache_states}")
    print(f"  Messages Exchanged: {messages}")
    print("-" * 40)

    # Operation 1: P1 reads X
    print("1. P1 reads X")
    # P1 issues a BusRd on a read miss.
    messages += 1
    message_counts.append(1)
    cache_states['P1'] = 'E'
    print("  Action: Read miss. P1 issues BusRd. State -> Exclusive (E).")
    print(f"  New States: {cache_states}")
    print(f"  Total Messages: {messages}")
    print("-" * 40)

    # Operation 2: P2 reads X
    print("2. P2 reads X")
    # P2 issues a BusRd. P1 snoops, provides data, and changes its state to Shared.
    messages += 1
    message_counts.append(1)
    cache_states['P1'] = 'S'
    cache_states['P2'] = 'S'
    print("  Action: Read miss. P2 issues BusRd. P1 provides data. States -> Shared (S).")
    print(f"  New States: {cache_states}")
    print(f"  Total Messages: {messages}")
    print("-" * 40)

    # Operation 3: P1 writes X = 1
    print("3. P1 writes X = 1")
    # P1 issues a BusRdX/Invalidate to upgrade its permission. P2 invalidates its copy.
    messages += 1
    message_counts.append(1)
    cache_states['P1'] = 'M'
    cache_states['P2'] = 'I'
    print("  Action: Write hit (permission fault). P1 issues BusRdX. P1 -> Modified (M), P2 -> Invalid (I).")
    print(f"  New States: {cache_states}")
    print(f"  Total Messages: {messages}")
    print("-" * 40)

    # Operation 4: P3 reads X
    print("4. P3 reads X")
    # P3 issues a BusRd. P1 snoops, provides data, and changes its state to Shared.
    messages += 1
    message_counts.append(1)
    cache_states['P1'] = 'S'
    cache_states['P3'] = 'S'
    print("  Action: Read miss. P3 issues BusRd. P1 provides data. P1 -> Shared (S), P3 -> Shared (S).")
    print(f"  New States: {cache_states}")
    print(f"  Total Messages: {messages}")
    print("-" * 40)

    # Operation 5: P2 writes X = 2
    print("5. P2 writes X = 2")
    # P2 issues a BusRdX on a write miss. P1 and P3 invalidate their copies.
    messages += 1
    message_counts.append(1)
    cache_states['P1'] = 'I'
    cache_states['P3'] = 'I'
    cache_states['P2'] = 'M'
    print("  Action: Write miss. P2 issues BusRdX. P2 -> Modified (M), P1/P3 -> Invalid (I).")
    print(f"  New States: {cache_states}")
    print(f"  Total Messages: {messages}")
    print("-" * 40)

    # Operation 6: P4 reads X
    print("6. P4 reads X")
    # P4 issues a BusRd. P2 snoops, provides data, and changes its state to Shared.
    messages += 1
    message_counts.append(1)
    cache_states['P2'] = 'S'
    cache_states['P4'] = 'S'
    print("  Action: Read miss. P4 issues BusRd. P2 provides data. P2 -> Shared (S), P4 -> Shared (S).")
    print(f"  New States: {cache_states}")
    print(f"  Total Messages: {messages}")
    print("-" * 40)

    # Operation 7: P1 reads X
    print("7. P1 reads X")
    # P1 issues a BusRd on a read miss. A sharer provides data.
    messages += 1
    message_counts.append(1)
    cache_states['P1'] = 'S'
    print("  Action: Read miss. P1 issues BusRd. P1 -> Shared (S).")
    print(f"  New States: {cache_states}")
    print(f"  Total Messages: {messages}")
    print("-" * 40)

    # Final Calculation
    print("\n--- Final Calculation ---")
    equation = " + ".join(map(str, message_counts))
    print("The total number of cache coherence messages is the sum of messages from each step:")
    print(f"{equation} = {messages}")

solve_mesi_simulation()