def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a sequence of processor operations and counts
    the number of cache coherence messages.
    """
    # Initial state
    # Each processor's cache state for variable X. M, E, S, I
    states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    # Value of X in memory. The 'owner' processor is also tracked.
    value_owner = 'memory'
    value = 0
    message_count = 0

    print("--- Initial State ---")
    print(f"Cache States: {states}")
    print(f"Value of X is {value} (in {value_owner})")
    print(f"Messages: {message_count}\n")
    print("--- Operation Sequence ---")

    # 1. P1 reads X
    print("1. P1 reads X")
    message_count += 1
    # P1 has a read miss, issues BusRd. No other cache has it, so memory responds.
    # P1 state becomes Exclusive.
    states['P1'] = 'E'
    value_owner = 'P1'
    print(f"   Bus Event: P1 issues BusRd. Memory supplies data.")
    print(f"   Action: P1's state becomes E (Exclusive).")
    print(f"   Cache States: {states}")
    print(f"   Total Messages: {message_count}\n")

    # 2. P2 reads X
    print("2. P2 reads X")
    message_count += 1
    # P2 has a read miss, issues BusRd. P1 snoops and has the block in E state.
    # P1 provides the data (cache-to-cache) and transitions to S.
    # P2's state becomes Shared.
    states['P1'] = 'S'
    states['P2'] = 'S'
    print(f"   Bus Event: P2 issues BusRd. P1 snoops and supplies data.")
    print(f"   Action: P1 state E -> S. P2 state I -> S.")
    print(f"   Cache States: {states}")
    print(f"   Total Messages: {message_count}\n")

    # 3. P1 writes X = 1
    print("3. P1 writes X = 1")
    value = 1
    value_owner = 'P1'
    # P1 has the block in S state. To write, it must invalidate other copies.
    # P1 issues a BusUpgr (invalidate) signal.
    message_count += 1
    # P2 snoops and invalidates its copy.
    states['P2'] = 'I'
    # P1's state becomes Modified.
    states['P1'] = 'M'
    print(f"   Bus Event: P1 issues BusUpgr to invalidate other shared copies.")
    print(f"   Action: P1 state S -> M. P2 state S -> I. P1 sets X=1.")
    print(f"   Cache States: {states}")
    print(f"   Total Messages: {message_count}\n")

    # 4. P3 reads X
    print("4. P3 reads X")
    message_count += 1
    # P3 has a read miss, issues BusRd. P1 snoops and has the block in M state.
    # P1 writes the block to memory, provides data to P3, and transitions to S.
    states['P1'] = 'S'
    states['P3'] = 'S'
    print(f"   Bus Event: P3 issues BusRd. P1 snoops, writes back to memory, and supplies data.")
    print(f"   Action: P1 state M -> S. P3 state I -> S.")
    print(f"   Cache States: {states}")
    print(f"   Total Messages: {message_count}\n")

    # 5. P2 writes X = 2
    print("5. P2 writes X = 2")
    value = 2
    value_owner = 'P2'
    # P2 has the block in I state (write miss). Issues BusRdX.
    message_count += 1
    # P1 and P3 snoop and invalidate their copies.
    states['P1'] = 'I'
    states['P3'] = 'I'
    # P2 gets the block and its state becomes Modified.
    states['P2'] = 'M'
    print(f"   Bus Event: P2 issues BusRdX (Read with Intent to Modify).")
    print(f"   Action: P1 and P3 invalidate (S -> I). P2 state I -> M. P2 sets X=2.")
    print(f"   Cache States: {states}")
    print(f"   Total Messages: {message_count}\n")

    # 6. P4 reads X
    print("6. P4 reads X")
    message_count += 1
    # P4 has a read miss, issues BusRd. P2 snoops (state=M).
    # P2 writes back to memory, provides data to P4, and transitions to S.
    states['P2'] = 'S'
    states['P4'] = 'S'
    print(f"   Bus Event: P4 issues BusRd. P2 snoops, writes back to memory, and supplies data.")
    print(f"   Action: P2 state M -> S. P4 state I -> S.")
    print(f"   Cache States: {states}")
    print(f"   Total Messages: {message_count}\n")

    # 7. P1 reads X
    print("7. P1 reads X")
    # P1 has the block in I state (read miss). Issues BusRd.
    message_count += 1
    # P2 and P4 have the block in S state. One of them (e.g., P2) supplies the data.
    states['P1'] = 'S'
    print(f"   Bus Event: P1 issues BusRd. P2 (or P4) snoops and supplies data.")
    print(f"   Action: P1 state I -> S. All other states (S) remain.")
    print(f"   Cache States: {states}")
    print(f"   Total Messages: {message_count}\n")

    print("--- Final Calculation ---")
    final_equation = "1 + 1 + 1 + 1 + 1 + 1 + 1"
    print(f"Total messages exchanged: {final_equation} = {message_count}")
    return message_count

if __name__ == "__main__":
    total = solve_mesi_coherence()
    # The final output is wrapped for the platform.
    # The script itself already prints the step-by-step breakdown.
    print(f"<<<{total}>>>")