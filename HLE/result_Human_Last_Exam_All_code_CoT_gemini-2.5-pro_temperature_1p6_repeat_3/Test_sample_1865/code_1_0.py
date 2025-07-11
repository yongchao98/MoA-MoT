def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a sequence of operations and counts
    the number of cache coherence messages.
    """
    # Initialize cache states for variable X. I=Invalid, S=Shared, E=Exclusive, M=Modified.
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    message_count = 0
    # Store the number of messages generated at each step
    messages_per_step = []

    def print_status(step, operation, explanation):
        print(f"--- Step {step}: {operation} ---")
        print(explanation)
        print(f"New Cache States: {caches}")
        print(f"Messages so far: {sum(messages_per_step)}\n")

    print("Initial State:")
    print(f"Cache States: {caches}")
    print("Value of X in memory = 0\n")

    # 1. P1 reads X
    # P1 has a cache miss (I). It issues a BusRd.
    # No other cache has X, so memory supplies it. P1's state becomes Exclusive (E).
    message_count = 1
    messages_per_step.append(1)
    caches['P1'] = 'E'
    explanation = ("P1 has a read miss. It sends a 'BusRd' message.\n"
                   "Memory responds. P1's cache state for X becomes Exclusive (E).")
    print_status(1, "P1 reads X", explanation)

    # 2. P2 reads X
    # P2 has a cache miss (I). It issues a BusRd.
    # P1 snoops, sees the request, and holds the line in E. P1 supplies the data.
    # Both P1 and P2 now hold the line in the Shared (S) state.
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    explanation = ("P2 has a read miss. It sends a 'BusRd' message.\n"
                   "P1 snoops and supplies the data. P1's state changes E -> S. P2's state becomes S.")
    print_status(2, "P2 reads X", explanation)

    # 3. P1 writes X = 1
    # P1 writes to a line in state S. This is a write hit, but requires invalidating other copies.
    # P1 sends a BusUpgr (Invalidate) message.
    # P2 snoops, sees the message, and invalidates its copy (S -> I). P1's state becomes Modified (M).
    message_count += 1
    messages_per_step.append(1)
    caches['P2'] = 'I'
    caches['P1'] = 'M'
    explanation = ("P1 writes to a shared line. It sends a 'BusUpgr' (Invalidate) message.\n"
                   "P2 invalidates its copy (S -> I). P1's state becomes Modified (M).")
    print_status(3, "P1 writes X = 1", explanation)

    # 4. P3 reads X
    # P3 has a cache miss (I). It issues a BusRd.
    # P1 snoops and has the line in state M. P1 supplies the data (now value 1).
    # P1's state changes to Shared (S), and P3's state becomes Shared (S).
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    explanation = ("P3 has a read miss. It sends a 'BusRd' message.\n"
                   "P1 snoops and supplies the data. P1's state changes M -> S. P3's state becomes S.")
    print_status(4, "P3 reads X", explanation)

    # 5. P2 writes X = 2
    # P2 has a write miss (I). It sends a BusRdX (Read Exclusive) message.
    # P1 and P3 snoop, see the BusRdX, and invalidate their copies (S -> I).
    # Data is supplied, and P2's state becomes Modified (M).
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'I'
    caches['P3'] = 'I'
    caches['P2'] = 'M'
    explanation = ("P2 has a write miss. It sends a 'BusRdX' message.\n"
                   "P1 and P3 invalidate their copies (S -> I). P2's state becomes Modified (M).")
    print_status(5, "P2 writes X = 2", explanation)

    # 6. P4 reads X
    # P4 has a cache miss (I). It issues a BusRd.
    # P2 snoops, has the line in state M, and supplies the data (value 2).
    # P2's state changes to Shared (S), and P4's state becomes Shared (S).
    message_count += 1
    messages_per_step.append(1)
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    explanation = ("P4 has a read miss. It sends a 'BusRd' message.\n"
                   "P2 snoops and supplies the data. P2's state changes M -> S. P4's state becomes S.")
    print_status(6, "P4 reads X", explanation)

    # 7. P1 reads X
    # P1 has a cache miss (I). It issues a BusRd.
    # P2 and P4 have the line in state S. A cache-to-cache transfer occurs.
    # P1's state becomes Shared (S). P2 and P4 remain in S.
    message_count += 1
    messages_per_step.append(1)
    caches['P1'] = 'S'
    explanation = ("P1 has a read miss. It sends a 'BusRd' message.\n"
                   "P2 or P4 supplies the data. P1's state becomes Shared (S).")
    print_status(7, "P1 reads X", explanation)

    # Final calculation
    equation = " + ".join(map(str, messages_per_step))
    print("------------------------------------------")
    print("Final Calculation Summary:")
    print("Each operation that causes a bus transaction for reading data or invalidating other caches counts as one message.")
    print(f"The total number of coherence messages is the sum of messages from each step:")
    print(f"{equation} = {message_count}")

solve_mesi_coherence()
print("<<<7>>>")