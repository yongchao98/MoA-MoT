def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a given sequence of operations and
    calculates the total number of cache coherence messages.
    """
    # Cache states: I (Invalid), S (Shared), E (Exclusive), M (Modified)
    p_states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    messages = 0
    message_log = []

    def print_states(operation):
        print(f"After '{operation}':")
        # Sort by processor name for consistent output
        sorted_processors = sorted(p_states.keys())
        state_str = ", ".join([f"{p}: {p_states[p]}" for p in sorted_processors])
        print(f"  States: [{state_str}]")
        print(f"  Total Messages: {messages}\n")

    # Initial state
    print("Initial State:")
    print(f"  Value of X in Memory = 0")
    sorted_processors = sorted(p_states.keys())
    state_str = ", ".join([f"{p}: {p_states[p]}" for p in sorted_processors])
    print(f"  States: [{state_str}]")
    print(f"  Total Messages: {messages}\n")


    # 1. P1 reads X
    step_messages = 0
    # P1 has a read miss, issues a BusRd. Gets data from memory. State becomes Exclusive.
    step_messages += 1
    messages += step_messages
    p_states['P1'] = 'E'
    message_log.append(step_messages)
    print("1. P1 reads X")
    print("  - P1 has a Read Miss.")
    print("  - P1 issues a 'BusRd' message to fetch data from memory. (+1 message)")
    print("  - No other cache shares the data, so P1's state becomes Exclusive (E).")
    print_states("P1 reads X")


    # 2. P2 reads X
    step_messages = 0
    # P2 has a read miss, issues a BusRd.
    step_messages += 1
    # P1 snoops, sees the BusRd, and supplies the data. P1's state changes from E -> S.
    p_states['P1'] = 'S'
    p_states['P2'] = 'S'
    messages += step_messages
    message_log.append(step_messages)
    print("2. P2 reads X")
    print("  - P2 has a Read Miss.")
    print("  - P2 issues a 'BusRd' message. (+1 message)")
    print("  - P1 snoops the bus, provides the data, and changes its state from E to Shared (S).")
    print("  - P2 receives the data and sets its state to Shared (S).")
    print_states("P2 reads X")


    # 3. P1 writes X = 1
    step_messages = 0
    # P1 has the block in state S, but needs exclusive access to write.
    # P1 issues a BusUpgr (invalidation) message.
    step_messages += 1
    # P2 snoops, sees the BusUpgr, and invalidates its copy (S -> I).
    p_states['P2'] = 'I'
    p_states['P1'] = 'M' # P1's state becomes Modified.
    messages += step_messages
    message_log.append(step_messages)
    print("3. P1 writes X = 1")
    print("  - P1 has a Write Hit, but its line is Shared (S).")
    print("  - P1 issues a 'BusUpgr' (Invalidate) message to gain exclusive ownership. (+1 message)")
    print("  - P2 snoops and invalidates its copy (S -> I).")
    print("  - P1's state becomes Modified (M).")
    print_states("P1 writes X = 1")


    # 4. P3 reads X
    step_messages = 0
    # P3 has a read miss, issues a BusRd.
    step_messages += 1
    # P1 snoops, sees the BusRd, and its state is M.
    # P1 supplies the data and also writes it back to memory (a separate bus transaction).
    step_messages += 1
    p_states['P1'] = 'S' # P1's state becomes S.
    p_states['P3'] = 'S' # P3 gets the data and becomes S.
    messages += step_messages
    message_log.append(step_messages)
    print("4. P3 reads X")
    print("  - P3 has a Read Miss and issues a 'BusRd' message. (+1 message)")
    print("  - P1 snoops the bus and sees its copy is Modified (M).")
    print("  - P1 provides the data to P3 and also issues a 'Write-Back' to memory. (+1 message)")
    print("  - P1's state changes to Shared (S), and P3's state becomes Shared (S).")
    print_states("P3 reads X")


    # 5. P2 writes X = 2
    step_messages = 0
    # P2 has a write miss (state I), issues BusRdX to get data and ownership.
    step_messages += 1
    # P1 and P3 snoop, see BusRdX, and invalidate their copies (S -> I).
    p_states['P1'] = 'I'
    p_states['P3'] = 'I'
    p_states['P2'] = 'M' # P2 gets data, writes, becomes M.
    messages += step_messages
    message_log.append(step_messages)
    print("5. P2 writes X = 2")
    print("  - P2 has a Write Miss.")
    print("  - P2 issues a 'BusRdX' (Read Exclusive) message to get data and invalidate others. (+1 message)")
    print("  - P1 and P3 snoop and invalidate their copies (S -> I).")
    print("  - P2 receives data, performs the write, and its state becomes Modified (M).")
    print_states("P2 writes X = 2")


    # 6. P4 reads X
    step_messages = 0
    # P4 has a read miss, issues BusRd.
    step_messages += 1
    # P2 snoops (state M), supplies data, and writes back to memory.
    step_messages += 1
    p_states['P2'] = 'S'
    p_states['P4'] = 'S'
    messages += step_messages
    message_log.append(step_messages)
    print("6. P4 reads X")
    print("  - P4 has a Read Miss and issues a 'BusRd' message. (+1 message)")
    print("  - P2 snoops, sees its copy is Modified (M), and intervenes.")
    print("  - P2 provides data to P4 and issues a 'Write-Back' to memory. (+1 message)")
    print("  - P2's state becomes Shared (S), and P4's state becomes Shared (S).")
    print_states("P4 reads X")


    # 7. P1 reads X
    step_messages = 0
    # P1 has a read miss, issues BusRd.
    step_messages += 1
    # P2 and P4 are Shared. One supplies the data (cache-to-cache).
    p_states['P1'] = 'S'
    messages += step_messages
    message_log.append(step_messages)
    print("7. P1 reads X")
    print("  - P1 has a Read Miss.")
    print("  - P1 issues a 'BusRd' message. (+1 message)")
    print("  - P2 and P4 are Shared. One of them provides the data via cache-to-cache transfer.")
    print("  - P1's state becomes Shared (S).")
    print_states("P1 reads X")
    
    # Final Calculation
    print("---------------------------------------------")
    print("Final Calculation")
    print("---------------------------------------------")
    print("The number of messages for each step are:")
    for i, count in enumerate(message_log):
        print(f"Step {i+1}: {count} message(s)")
        
    equation = " + ".join(map(str, message_log))
    print(f"\nTotal Messages = {equation} = {messages}")


if __name__ == "__main__":
    solve_mesi_coherence()
    print("\n<<<9>>>")