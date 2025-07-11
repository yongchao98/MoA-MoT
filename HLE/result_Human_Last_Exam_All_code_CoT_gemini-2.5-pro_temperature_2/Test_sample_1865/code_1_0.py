def simulate_mesi_protocol():
    """
    Simulates the MESI protocol for a given sequence of operations
    and counts the cache coherence messages.
    """
    # Cache states for variable X for each processor [P1, P2, P3, P4]
    # States: M (Modified), E (Exclusive), S (Shared), I (Invalid)
    states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    messages = 0
    
    print("Initial State:")
    print(f"  Processor States: {states}")
    print(f"  Coherence Messages: {messages}\n" + "-"*60)

    # Operation 1: P1 reads X
    # P1 has a read miss, issues a BusRd. No other cache has X, so P1 takes it in Exclusive state.
    messages += 1
    states['P1'] = 'E'
    print("Step 1: P1 reads X")
    print("  - P1's cache is 'Invalid', causing a read miss.")
    print("  - P1 issues a 'BusRd' message to read from memory.")
    print(f"  - Message count is now {messages}.")
    print(f"  - Since P1 is the only cache with X, its state becomes 'Exclusive' (E).")
    print(f"  Processor States: {states}\n" + "-"*60)

    # Operation 2: P2 reads X
    # P2 has a read miss, issues a BusRd. P1 snoops, sees the read, provides the data, and changes its state to Shared.
    # P2 takes the data in Shared state.
    messages += 1
    states['P1'] = 'S'
    states['P2'] = 'S'
    print("Step 2: P2 reads X")
    print("  - P2's cache is 'Invalid', causing a read miss.")
    print("  - P2 issues a 'BusRd' message.")
    print(f"  - Message count is now {messages}.")
    print("  - P1 snoops the bus, provides the data to P2, and changes its state from 'E' to 'S'.")
    print("  - P2 receives the data and sets its state to 'Shared' (S).")
    print(f"  Processor States: {states}\n" + "-"*60)

    # Operation 3: P1 writes X = 1
    # P1 has X in Shared state. A write requires an upgrade to Modified.
    # P1 issues an Invalidate message to invalidate all other shared copies (P2).
    messages += 1
    states['P1'] = 'M'
    states['P2'] = 'I'
    print("Step 3: P1 writes X = 1")
    print("  - P1 holds the block in 'Shared' (S), but a write needs exclusive access.")
    print("  - P1 issues an 'Invalidate' message to invalidate other copies.")
    print(f"  - Message count is now {messages}.")
    print("  - P2 receives the invalidate message and changes its state from 'S' to 'I'.")
    print("  - P1's state becomes 'Modified' (M).")
    print(f"  Processor States: {states}\n" + "-"*60)
    
    # Operation 4: P3 reads X
    # P3 has a read miss, issues BusRd. P1 snoops, sees it has the block in Modified state.
    # P1 flushes the data to the bus (for P3 and memory) and changes its state to Shared. P3 becomes Shared.
    messages += 1
    states['P1'] = 'S'
    states['P3'] = 'S'
    print("Step 4: P3 reads X")
    print("  - P3's cache is 'Invalid', causing a read miss.")
    print("  - P3 issues a 'BusRd' message.")
    print(f"  - Message count is now {messages}.")
    print("  - P1 snoops the bus, provides the modified data, and changes its state from 'M' to 'S'.")
    print("  - P3 receives the data and sets its state to 'Shared' (S).")
    print(f"  Processor States: {states}\n" + "-"*60)
    
    # Operation 5: P2 writes X = 2
    # P2 has X in Invalid state, a write miss.
    # P2 issues a BusRdX (Read for Ownership) message to get data and invalidate others.
    # P1 and P3 invalidate their copies. P2 becomes Modified.
    messages += 1
    states['P1'] = 'I'
    states['P3'] = 'I'
    states['P2'] = 'M'
    print("Step 5: P2 writes X = 2")
    print("  - P2's cache is 'Invalid', causing a write miss.")
    print("  - P2 issues a 'BusRdX' (Read for Ownership) message.")
    print(f"  - Message count is now {messages}.")
    print("  - P1 and P3 snoop the BusRdX, invalidate their copies ('S' -> 'I').")
    print("  - P2 receives the data and its state becomes 'Modified' (M) after writing.")
    print(f"  Processor States: {states}\n" + "-"*60)
    
    # Operation 6: P4 reads X
    # P4 has a read miss, issues BusRd. P2 snoops, has it in Modified.
    # P2 flushes the data and changes to Shared. P4 becomes Shared.
    messages += 1
    states['P2'] = 'S'
    states['P4'] = 'S'
    print("Step 6: P4 reads X")
    print("  - P4's cache is 'Invalid', causing a read miss.")
    print("  - P4 issues a 'BusRd' message.")
    print(f"  - Message count is now {messages}.")
    print("  - P2 snoops the bus, provides the modified data, and changes its state from 'M' to 'S'.")
    print("  - P4 receives the data and sets its state to 'Shared' (S).")
    print(f"  Processor States: {states}\n" + "-"*60)

    # Operation 7: P1 reads X
    # P1 has a read miss, issues BusRd. P2 and P4 have the data in Shared state.
    # Data is provided by one of the caches (or memory, since it's clean). No state changes for P2/P4.
    # P1's state becomes Shared.
    messages += 1
    states['P1'] = 'S'
    print("Step 7: P1 reads X")
    print("  - P1's cache is 'Invalid', causing a read miss.")
    print("  - P1 issues a 'BusRd' message.")
    print(f"  - Message count is now {messages}.")
    print("  - P2 and P4 have the data in 'Shared' (S) state. One of them supplies the data.")
    print("  - P1 receives the data and sets its state to 'Shared' (S).")
    print(f"  Processor States: {states}\n" + "-"*60)

    print("Final Result:")
    print(f"The total number of cache coherence messages exchanged is {messages}.")
    
if __name__ == "__main__":
    simulate_mesi_protocol()