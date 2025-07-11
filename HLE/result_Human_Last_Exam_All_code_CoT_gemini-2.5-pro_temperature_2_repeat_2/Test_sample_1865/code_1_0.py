def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a sequence of operations and counts
    the cache coherence messages exchanged.
    """
    # Initial system state
    # Each cache holds the state ('I', 'S', 'E', 'M') and the value of X
    caches = {f'P{i}': {'state': 'I', 'value': None} for i in range(1, 5)}
    memory_value = 0
    
    # List to track the count of messages at each step for the final equation
    message_counts = []
    
    print("### MESI Protocol Simulation ###\n")
    print("--- Initial State ---")
    print(f"Caches: {[c['state'] for c in caches.values()]}, Memory X: {memory_value}\n")

    # 1. P1 reads X
    # P1 has a Read Miss. It sends a 'Read Miss' message.
    # Since no other cache has X, memory supplies it. P1's state becomes Exclusive.
    message_counts.append("1")
    caches['P1']['state'] = 'E'
    caches['P1']['value'] = memory_value
    print("1. P1 reads X")
    print("   - P1 issues a 'Read Miss'. (1 message)")
    print("   - P1 cache state becomes 'Exclusive (E)'.")
    print(f"   - Current States: {[c['state'] for c in caches.values()]}, Message Count: {len(message_counts)}\n")

    # 2. P2 reads X
    # P2 has a Read Miss. It sends a 'Read Miss' message.
    # P1 snoops the bus, sees the request, and supplies the data (cache-to-cache).
    # P1's state changes from Exclusive to Shared. P2's state becomes Shared.
    message_counts.append("1")
    caches['P1']['state'] = 'S'
    caches['P2']['state'] = 'S'
    caches['P2']['value'] = caches['P1']['value']
    print("2. P2 reads X")
    print("   - P2 issues a 'Read Miss'. (1 message)")
    print("   - P1 supplies data. P1 state -> 'Shared (S)'.")
    print("   - P2 cache state becomes 'Shared (S)'.")
    print(f"   - Current States: {[c['state'] for c in caches.values()]}, Message Count: {len(message_counts)}\n")

    # 3. P1 writes X = 1
    # P1 is in Shared state. To write, it must upgrade to Modified.
    # It sends an 'Upgrade' (or 'Read for Ownership') message to invalidate other copies.
    # P2 receives this message and invalidates its copy (S -> I).
    # P1's state becomes Modified.
    message_counts.append("1")
    caches['P1']['state'] = 'M'
    caches['P1']['value'] = 1
    caches['P2']['state'] = 'I'
    caches['P2']['value'] = None
    print("3. P1 writes X = 1")
    print("   - P1 issues an 'Upgrade' message to invalidate sharers. (1 message)")
    print("   - P2 invalidates its cache. P1 state -> 'Modified (M)'.")
    print(f"   - Current States: {[c['state'] for c in caches.values()]}, Message Count: {len(message_counts)}\n")

    # 4. P3 reads X
    # P3 has a Read Miss. It sends a 'Read Miss' message.
    # P1 has the data in Modified state. It writes the data back to memory and supplies it to P3.
    # P1's state changes from Modified to Shared. P3's state becomes Shared.
    message_counts.append("1")
    memory_value = caches['P1']['value'] # P1 writes back to memory
    caches['P1']['state'] = 'S'
    caches['P3']['state'] = 'S'
    caches['P3']['value'] = memory_value
    print("4. P3 reads X")
    print("   - P3 issues a 'Read Miss'. (1 message)")
    print("   - P1 flushes data to memory and P3. P1 state -> 'Shared (S)'.")
    print("   - P3 cache state becomes 'Shared (S)'.")
    print(f"   - Current States: {[c['state'] for c in caches.values()]}, Message Count: {len(message_counts)}\n")

    # 5. P2 writes X = 2
    # P2 has a Write Miss (its state is Invalid).
    # It sends a 'Read for Ownership' message.
    # P1 and P3 see the message and invalidate their copies (S -> I).
    # P2 gets the data, performs the write, and its state becomes Modified.
    message_counts.append("1")
    caches['P1']['state'] = 'I'
    caches['P1']['value'] = None
    caches['P3']['state'] = 'I'
    caches['P3']['value'] = None
    caches['P2']['state'] = 'M'
    caches['P2']['value'] = 2
    print("5. P2 writes X = 2")
    print("   - P2 issues a 'Read for Ownership' message. (1 message)")
    print("   - P1 and P3 invalidate. P2 state -> 'Modified (M)'.")
    print(f"   - Current States: {[c['state'] for c in caches.values()]}, Message Count: {len(message_counts)}\n")

    # 6. P4 reads X
    # P4 has a Read Miss. It sends a 'Read Miss' message.
    # P2 has the data in Modified state. It writes the data back to memory and supplies it to P4.
    # P2's state changes from Modified to Shared. P4's state becomes Shared.
    message_counts.append("1")
    memory_value = caches['P2']['value'] # P2 writes back to memory
    caches['P2']['state'] = 'S'
    caches['P4']['state'] = 'S'
    caches['P4']['value'] = memory_value
    print("6. P4 reads X")
    print("   - P4 issues a 'Read Miss'. (1 message)")
    print("   - P2 flushes data to memory and P4. P2 state -> 'Shared (S)'.")
    print("   - P4 cache state becomes 'Shared (S)'.")
    print(f"   - Current States: {[c['state'] for c in caches.values()]}, Message Count: {len(message_counts)}\n")

    # 7. P1 reads X
    # P1 has a Read Miss. It sends a 'Read Miss' message.
    # P2 and P4 are sharers. One of them supplies the data to P1.
    # P1's state becomes Shared. States of P2 and P4 are unchanged.
    message_counts.append("1")
    caches['P1']['state'] = 'S'
    caches['P1']['value'] = caches['P2']['value']
    print("7. P1 reads X")
    print("   - P1 issues a 'Read Miss'. (1 message)")
    print("   - A sharer (P2 or P4) supplies data. P1 state -> 'Shared (S)'.")
    print(f"   - Current States: {[c['state'] for c in caches.values()]}, Message Count: {len(message_counts)}\n")

    # Final count
    total_messages = len(message_counts)
    equation = " + ".join(message_counts)
    print("--- Final Calculation ---")
    print(f"The total number of messages is the sum of messages from each step.")
    print(f"Calculation: {equation} = {total_messages}")
    print(f"\nTotal Cache Coherence Messages: {total_messages}")
    print(f"<<<{total_messages}>>>")

# Run the simulation
solve_mesi_coherence()