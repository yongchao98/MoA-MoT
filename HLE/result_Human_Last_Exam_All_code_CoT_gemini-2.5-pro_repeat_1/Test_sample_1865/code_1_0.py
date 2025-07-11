def solve_mesi_problem():
    """
    Simulates the MESI protocol for a given sequence of operations and
    counts the number of cache coherence messages.
    """
    # Initial state of the system
    caches = {
        'P1': {'state': 'I', 'value': None},
        'P2': {'state': 'I', 'value': None},
        'P3': {'state': 'I', 'value': None},
        'P4': {'state': 'I', 'value': None}
    }
    memory_value = 0
    message_count = 0
    message_components = []

    def print_state(operation_description):
        print(f"--- {operation_description} ---")
        print(f"Cache States: P1:{caches['P1']['state']}, P2:{caches['P2']['state']}, P3:{caches['P3']['state']}, P4:{caches['P4']['state']}")
        print(f"Memory[X] = {memory_value}")
        print(f"Current Message Count: {message_count}\n")

    print("Initial State:")
    print(f"X is not cached. Memory[X] = {memory_value}\n")

    # --- Step 1: P1 reads X ---
    # P1 has a read miss (state is Invalid).
    # P1 issues a BusRd message to the bus.
    # Data is supplied by main memory. Since no other cache has the block, P1's state becomes Exclusive (E).
    message_count += 1
    message_components.append(1)
    caches['P1']['state'] = 'E'
    caches['P1']['value'] = memory_value
    print_state("Step 1: P1 reads X")

    # --- Step 2: P2 reads X ---
    # P2 has a read miss (state is Invalid).
    # P2 issues a BusRd message.
    # P1 snoops the bus, sees the request for a block it holds in Exclusive (E) state.
    # P1 supplies the data to P2 and changes its state to Shared (S).
    # P2 receives the data and sets its state to Shared (S).
    message_count += 1
    message_components.append(1)
    caches['P1']['state'] = 'S'
    caches['P2']['state'] = 'S'
    caches['P2']['value'] = caches['P1']['value']
    print_state("Step 2: P2 reads X")

    # --- Step 3: P1 writes X = 1 ---
    # P1 wants to write, but its state is Shared (S). It must gain exclusive ownership.
    # P1 issues a BusRdX (or BusUpgr) message to invalidate other copies.
    # P2 snoops the bus, invalidates its copy (S -> I).
    # P1's state changes to Modified (M) and it writes the new value.
    message_count += 1
    message_components.append(1)
    caches['P2']['state'] = 'I'
    caches['P2']['value'] = None
    caches['P1']['state'] = 'M'
    caches['P1']['value'] = 1
    print_state("Step 3: P1 writes X = 1")

    # --- Step 4: P3 reads X ---
    # P3 has a read miss (state is Invalid).
    # P3 issues a BusRd message.
    # P1 snoops the bus, sees the request for a block it holds in Modified (M) state.
    # P1 supplies the data (1) to P3 and writes the value back to main memory.
    # P1's state changes from Modified (M) to Shared (S).
    # P3 receives the data and sets its state to Shared (S).
    message_count += 1
    message_components.append(1)
    memory_value = caches['P1']['value']
    caches['P1']['state'] = 'S'
    caches['P3']['state'] = 'S'
    caches['P3']['value'] = memory_value
    print_state("Step 4: P3 reads X")

    # --- Step 5: P2 writes X = 2 ---
    # P2 has a write miss (state is Invalid).
    # P2 issues a BusRdX message to get the data and invalidate others.
    # P1 and P3 snoop the bus, see the request, and invalidate their copies (S -> I).
    # Data is supplied (e.g., from P1's cache), and P2's state becomes Modified (M). P2 writes the new value.
    message_count += 1
    message_components.append(1)
    caches['P1']['state'] = 'I'
    caches['P1']['value'] = None
    caches['P3']['state'] = 'I'
    caches['P3']['value'] = None
    caches['P2']['state'] = 'M'
    caches['P2']['value'] = 2
    print_state("Step 5: P2 writes X = 2")

    # --- Step 6: P4 reads X ---
    # P4 has a read miss (state is Invalid).
    # P4 issues a BusRd message.
    # P2 snoops the bus, sees the request for a block it holds in Modified (M) state.
    # P2 supplies the data (2) to P4 and writes the value back to main memory.
    # P2's state changes from Modified (M) to Shared (S).
    # P4 receives the data and sets its state to Shared (S).
    message_count += 1
    message_components.append(1)
    memory_value = caches['P2']['value']
    caches['P2']['state'] = 'S'
    caches['P4']['state'] = 'S'
    caches['P4']['value'] = memory_value
    print_state("Step 6: P4 reads X")

    # --- Step 7: P1 reads X ---
    # P1 has a read miss (state is Invalid).
    # P1 issues a BusRd message.
    # P2 and P4 have the data in Shared (S) state. One of them (e.g., P2) supplies the data.
    # P1 receives the data and sets its state to Shared (S). All sharers remain in state S.
    message_count += 1
    message_components.append(1)
    caches['P1']['state'] = 'S'
    caches['P1']['value'] = caches['P2']['value']
    print_state("Step 7: P1 reads X")

    # --- Final Calculation ---
    equation = " + ".join(map(str, message_components))
    print("--- Final Calculation ---")
    print(f"The total number of cache coherence messages is the sum of messages from each step:")
    print(f"{equation} = {message_count}")


if __name__ == '__main__':
    solve_mesi_problem()
    print("\n<<<7>>>")