def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a given sequence of operations and
    counts the number of cache coherence messages.
    """
    # Initialize cache states for variable X. I=Invalid, E=Exclusive, S=Shared, M=Modified.
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    total_messages = 0
    message_log = []

    print("Initial State:")
    print(f"  - Caches: {caches}")
    print("-" * 50)

    # --- Step 1: P1 reads X ---
    print("Step 1: P1 reads X")
    print(f"  - Initial States: {caches}")
    # P1 has a read miss. It sends a BusRd message to memory.
    # No other cache has the data, so P1 loads it in the Exclusive state.
    total_messages += 1
    message_log.append('1')
    caches['P1'] = 'E'
    print("  - Action: P1 has a read miss, sends BusRd message.")
    print("  - State Change: P1 loads X from memory and sets its state to Exclusive (E).")
    print(f"  - Final States: {caches}")
    print(f"  - Messages for this step: 1. Total Messages: {total_messages}")
    print("-" * 50)

    # --- Step 2: P2 reads X ---
    print("Step 2: P2 reads X")
    print(f"  - Initial States: {caches}")
    # P2 has a read miss. It sends a BusRd message.
    # P1 snoops the bus, sees it has the data (in state E), provides the data, and changes its state to Shared.
    # P2 loads the data and sets its state to Shared.
    total_messages += 1
    message_log.append('1')
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print("  - Action: P2 has a read miss, sends BusRd message.")
    print("  - State Change: P1 provides the data and changes to Shared (S). P2 loads data and becomes Shared (S).")
    print(f"  - Final States: {caches}")
    print(f"  - Messages for this step: 1. Total Messages: {total_messages}")
    print("-" * 50)

    # --- Step 3: P1 writes X = 1 ---
    print("Step 3: P1 writes X = 1")
    print(f"  - Initial States: {caches}")
    # P1 wants to write but is in Shared state. It must invalidate other copies.
    # P1 sends a BusUpgr (or Invalidate) message on the bus.
    # P2 snoops the bus and invalidates its copy. P1's state becomes Modified.
    total_messages += 1
    message_log.append('1')
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print("  - Action: P1 has a write hit, but needs to invalidate others. Sends BusUpgr message.")
    print("  - State Change: P1 updates its state to Modified (M). P2 invalidates its copy to Invalid (I).")
    print(f"  - Final States: {caches}")
    print(f"  - Messages for this step: 1. Total Messages: {total_messages}")
    print("-" * 50)

    # --- Step 4: P3 reads X ---
    print("Step 4: P3 reads X")
    print(f"  - Initial States: {caches}")
    # P3 has a read miss. It sends a BusRd message.
    # P1 snoops the bus, sees it has the modified data, provides it, and writes it back to memory.
    # P1's state changes to Shared. P3 loads the data and becomes Shared.
    total_messages += 1
    message_log.append('1')
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print("  - Action: P3 has a read miss, sends BusRd message.")
    print("  - State Change: P1 provides the modified data, writes back to memory, and becomes Shared (S). P3 becomes Shared (S).")
    print(f"  - Final States: {caches}")
    print(f"  - Messages for this step: 1. Total Messages: {total_messages}")
    print("-" * 50)

    # --- Step 5: P2 writes X = 2 ---
    print("Step 5: P2 writes X = 2")
    print(f"  - Initial States: {caches}")
    # P2 has a write miss (its state is Invalid).
    # P2 sends a BusRdX (Read Exclusive) message to get the data and invalidate others.
    # P1 and P3 snoop, invalidate their copies. P2 gets the data, writes its value, and becomes Modified.
    total_messages += 1
    message_log.append('1')
    caches['P2'] = 'M'
    caches['P1'] = 'I'
    caches['P3'] = 'I'
    print("  - Action: P2 has a write miss, sends BusRdX message to acquire exclusive ownership.")
    print("  - State Change: P2 becomes Modified (M). P1 and P3 invalidate their copies to Invalid (I).")
    print(f"  - Final States: {caches}")
    print(f"  - Messages for this step: 1. Total Messages: {total_messages}")
    print("-" * 50)

    # --- Step 6: P4 reads X ---
    print("Step 6: P4 reads X")
    print(f"  - Initial States: {caches}")
    # P4 has a read miss. It sends a BusRd message.
    # P2 snoops the bus, provides the modified data, writes back to memory.
    # P2's state becomes Shared. P4 loads data and becomes Shared.
    total_messages += 1
    message_log.append('1')
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print("  - Action: P4 has a read miss, sends BusRd message.")
    print("  - State Change: P2 provides data, writes back to memory, and becomes Shared (S). P4 becomes Shared (S).")
    print(f"  - Final States: {caches}")
    print(f"  - Messages for this step: 1. Total Messages: {total_messages}")
    print("-" * 50)
    
    # --- Step 7: P1 reads X ---
    print("Step 7: P1 reads X")
    print(f"  - Initial States: {caches}")
    # P1 has a read miss (state is Invalid).
    # P1 sends a BusRd message.
    # P2 and P4 are Shared. One of them provides the data. States remain Shared. P1 becomes Shared.
    total_messages += 1
    message_log.append('1')
    caches['P1'] = 'S'
    print("  - Action: P1 has a read miss, sends BusRd message.")
    print("  - State Change: P2 (or P4) provides the data. P1 loads data and becomes Shared (S).")
    print(f"  - Final States: {caches}")
    print(f"  - Messages for this step: 1. Total Messages: {total_messages}")
    print("-" * 50)
    
    print("\nSummary of Coherence Messages:")
    equation_str = " + ".join(message_log)
    print(f"Each of the 7 operations required one bus message to maintain coherence.")
    print(f"The total number of messages is the sum of messages from each step.")
    print(f"\nFinal Equation: {equation_str} = {total_messages}")


if __name__ == '__main__':
    solve_mesi_coherence()
    print("\n<<<7>>>")