def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a given sequence of operations
    and calculates the total number of coherence messages.
    """
    # Initial state
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    memory_value = 0
    message_count = 0
    total_messages = 0

    print("--- MESI Protocol Simulation ---")
    print(f"Initial State: Caches={caches}, Memory_X={memory_value}, Messages={total_messages}\n")

    # --- Operation 1: P1 reads X ---
    message_count = 1
    total_messages += message_count
    caches['P1'] = 'E'
    print("1. P1 reads X")
    print("   - P1 has a cache miss (I -> E).")
    print("   - P1 sends a 'Read Miss' message to the bus.")
    print(f"   - Messages: {message_count}")
    print(f"   - State: Caches={caches}, Messages so far={total_messages}\n")
    # P1 + 1 = 1
    op1_msg = 1
    
    # --- Operation 2: P2 reads X ---
    message_count = 1
    total_messages += message_count
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print("2. P2 reads X")
    print("   - P2 has a cache miss.")
    print("   - P2 sends a 'Read Miss' message. P1 snoops and provides the data.")
    print("   - P1's state changes E -> S. P2's state becomes S.")
    print(f"   - Messages: {message_count}")
    print(f"   - State: Caches={caches}, Messages so far={total_messages}\n")
    # P2 + 1 = 2
    op2_msg = 1
    
    # --- Operation 3: P1 writes X = 1 ---
    message_count = 1
    total_messages += message_count
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    memory_value = "stale(0)"
    print("3. P1 writes X = 1")
    print("   - P1 has the block in state S. It must upgrade to write.")
    print("   - P1 sends an 'Invalidate' message on the bus.")
    print("   - P1's state becomes M. P2 invalidates its copy (S -> I).")
    print(f"   - Messages: {message_count}")
    print(f"   - State: Caches={caches}, Messages so far={total_messages}\n")
    # P1 + 1 = 3
    op3_msg = 1
    
    # --- Operation 4: P3 reads X ---
    message_count = 1
    total_messages += message_count
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    memory_value = 1
    print("4. P3 reads X")
    print("   - P3 has a cache miss and sends a 'Read Miss' message.")
    print("   - P1 (in state M) provides the data, writes it to memory, and changes its state to S.")
    print("   - P3's state becomes S.")
    print(f"   - Messages: {message_count}")
    print(f"   - State: Caches={caches}, Messages so far={total_messages}\n")
    # P3 + 1 = 4
    op4_msg = 1
    
    # --- Operation 5: P2 writes X = 2 ---
    message_count = 1
    total_messages += message_count
    caches['P2'] = 'M'
    caches['P1'] = 'I'
    caches['P3'] = 'I'
    memory_value = "stale(1)"
    print("5. P2 writes X = 2")
    print("   - P2 has a write miss (state I) and sends a 'Read For Ownership' message.")
    print("   - P1 and P3 invalidate their copies (S -> I). P2's state becomes M.")
    print(f"   - Messages: {message_count}")
    print(f"   - State: Caches={caches}, Messages so far={total_messages}\n")
    # P2 + 1 = 5
    op5_msg = 1

    # --- Operation 6: P4 reads X ---
    message_count = 1
    total_messages += message_count
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    memory_value = 2
    print("6. P4 reads X")
    print("   - P4 has a cache miss and sends a 'Read Miss' message.")
    print("   - P2 (in state M) provides the data, writes it to memory, and changes its state to S.")
    print("   - P4's state becomes S.")
    print(f"   - Messages: {message_count}")
    print(f"   - State: Caches={caches}, Messages so far={total_messages}\n")
    # P4 + 1 = 6
    op6_msg = 1

    # --- Operation 7: P1 reads X ---
    message_count = 1
    total_messages += message_count
    caches['P1'] = 'S'
    print("7. P1 reads X")
    print("   - P1 has a cache miss and sends a 'Read Miss' message.")
    print("   - P2 or P4 (in state S) provides the data.")
    print("   - P1's state becomes S.")
    print(f"   - Messages: {message_count}")
    print(f"   - State: Caches={caches}, Messages so far={total_messages}\n")
    # P1 + 1 = 7
    op7_msg = 1
    
    print("--- Calculation Summary ---")
    print(f"The total number of cache coherence messages is the sum of messages from each step:")
    print(f"{op1_msg} + {op2_msg} + {op3_msg} + {op4_msg} + {op5_msg} + {op6_msg} + {op7_msg} = {total_messages}")

solve_mesi_coherence()
<<<7>>>