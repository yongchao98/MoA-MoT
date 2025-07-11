def solve_mesi_messages():
    """
    Simulates the MESI protocol for a sequence of operations and counts the
    coherence messages.

    A coherence message is defined as a bus transaction initiated by a processor
    to maintain cache coherence (Read Miss, Invalidate, Read For Ownership).
    """

    # Initial states
    caches = {
        'P1': {'state': 'I', 'value': None},
        'P2': {'state': 'I', 'value': None},
        'P3': {'state': 'I', 'value': None},
        'P4': {'state': 'I', 'value': None},
    }
    memory_value = 0
    message_counts = []

    print("--- Initial State ---")
    print(f"Caches: P1={caches['P1']['state']}, P2={caches['P2']['state']}, P3={caches['P3']['state']}, P4={caches['P4']['state']}")
    print(f"Memory X = {memory_value}\n")

    # Operation 1: P1 reads X
    print("1. P1 reads X")
    message_counts.append(1)
    # P1 has a Read Miss, broadcasts Read Miss on the bus (1 message)
    # Memory responds. P1 is the only cache, so state becomes Exclusive.
    caches['P1']['state'] = 'E'
    caches['P1']['value'] = memory_value
    print("   -> P1 Miss. Issues Read Miss (1 message). State -> E.")
    print(f"   Caches: P1={caches['P1']['state']}, P2={caches['P2']['state']}, P3={caches['P3']['state']}, P4={caches['P4']['state']}\n")

    # Operation 2: P2 reads X
    print("2. P2 reads X")
    message_counts.append(1)
    # P2 has a Read Miss, broadcasts Read Miss (1 message)
    # P1 snoops, has X in state E. P1 provides data and transitions to Shared.
    caches['P1']['state'] = 'S'
    caches['P2']['state'] = 'S'
    caches['P2']['value'] = caches['P1']['value']
    print("   -> P2 Miss. Issues Read Miss (1 message). P1 supplies data, P1: E->S, P2 -> S.")
    print(f"   Caches: P1={caches['P1']['state']}, P2={caches['P2']['state']}, P3={caches['P3']['state']}, P4={caches['P4']['state']}\n")

    # Operation 3: P1 writes X = 1
    print("3. P1 writes X = 1")
    message_counts.append(1)
    # P1 is in Shared state, needs to upgrade. Broadcasts Invalidate (1 message).
    # P2 snoops and invalidates its copy.
    caches['P2']['state'] = 'I'
    # P1 writes and transitions to Modified.
    caches['P1']['state'] = 'M'
    caches['P1']['value'] = 1
    print("   -> P1 Write Hit (in S). Issues Invalidate (1 message). P2: S->I, P1: S->M.")
    print(f"   Caches: P1={caches['P1']['state']}, P2={caches['P2']['state']}, P3={caches['P3']['state']}, P4={caches['P4']['state']}\n")

    # Operation 4: P3 reads X
    print("4. P3 reads X")
    message_counts.append(1)
    # P3 has a Read Miss, broadcasts Read Miss (1 message)
    # P1 snoops (state M), provides data, writes back to memory. P1 transitions to Shared.
    memory_value = caches['P1']['value']
    caches['P1']['state'] = 'S'
    caches['P3']['state'] = 'S'
    caches['P3']['value'] = caches['P1']['value']
    print("   -> P3 Miss. Issues Read Miss (1 message). P1 supplies data, P1: M->S, P3 -> S.")
    print(f"   Caches: P1={caches['P1']['state']}, P2={caches['P2']['state']}, P3={caches['P3']['state']}, P4={caches['P4']['state']}\n")

    # Operation 5: P2 writes X = 2
    print("5. P2 writes X = 2")
    message_counts.append(1)
    # P2 has a Write Miss, broadcasts RFO (Read For Ownership) (1 message).
    # P1 and P3 snoop, invalidate their copies (S->I).
    caches['P1']['state'] = 'I'
    caches['P3']['state'] = 'I'
    # P2 gets data, writes, and transitions to Modified.
    caches['P2']['state'] = 'M'
    caches['P2']['value'] = 2
    print("   -> P2 Write Miss. Issues RFO (1 message). P1/P3 invalidate. P2 -> M.")
    print(f"   Caches: P1={caches['P1']['state']}, P2={caches['P2']['state']}, P3={caches['P3']['state']}, P4={caches['P4']['state']}\n")

    # Operation 6: P4 reads X
    print("6. P4 reads X")
    message_counts.append(1)
    # P4 has a Read Miss, broadcasts Read Miss (1 message).
    # P2 snoops (state M), provides data, writes back. P2 transitions to Shared.
    memory_value = caches['P2']['value']
    caches['P2']['state'] = 'S'
    caches['P4']['state'] = 'S'
    caches['P4']['value'] = caches['P2']['value']
    print("   -> P4 Miss. Issues Read Miss (1 message). P2 supplies data, P2: M->S, P4 -> S.")
    print(f"   Caches: P1={caches['P1']['state']}, P2={caches['P2']['state']}, P3={caches['P3']['state']}, P4={caches['P4']['state']}\n")

    # Operation 7: P1 reads X
    print("7. P1 reads X")
    message_counts.append(1)
    # P1 has a Read Miss, broadcasts Read Miss (1 message).
    # P2/P4 are in state S and one provides data.
    caches['P1']['state'] = 'S'
    caches['P1']['value'] = caches['P2']['value']
    print("   -> P1 Miss. Issues Read Miss (1 message). P2/P4 supply data, P1 -> S.")
    print(f"   Caches: P1={caches['P1']['state']}, P2={caches['P2']['state']}, P3={caches['P3']['state']}, P4={caches['P4']['state']}\n")

    # Final calculation
    total_messages = sum(message_counts)
    equation = " + ".join(map(str, message_counts))
    print("--- Final Calculation ---")
    print(f"Total Messages = {equation} = {total_messages}")

solve_mesi_messages()
<<<7>>>