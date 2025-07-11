def solve_mesi_coherence():
    """
    Calculates the number of cache coherence messages for a given sequence of operations
    in a 4-processor system using the MESI protocol.
    """
    # Initial states
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    memory_value = 0
    message_count = 0
    step_messages = []

    def print_state():
        print(f"Current States: P1:{caches['P1']}, P2:{caches['P2']}, P3:{caches['P3']}, P4:{caches['P4']}")
        print("-" * 40)

    print("Initial State: Variable X is not in any cache, value in memory is 0.\n")
    print_state()

    # --- Step 1: P1 reads X ---
    print("1. P1 reads X")
    print("   - P1 has a cache miss.")
    print("   - P1 sends a 'BusRd' message to the bus. (1 message)")
    print("   - No other cache has X. Memory responds with the data. (1 message)")
    print("   - P1 loads X in the 'Exclusive' (E) state as it's the only copy.")
    messages_this_step = 2
    step_messages.append(messages_this_step)
    message_count += messages_this_step
    caches['P1'] = 'E'
    print(f"   Messages exchanged in this step: {messages_this_step}")
    print(f"   Total messages so far: {message_count}")
    print_state()

    # --- Step 2: P2 reads X ---
    print("2. P2 reads X")
    print("   - P2 has a cache miss.")
    print("   - P2 sends a 'BusRd' message to the bus. (1 message)")
    print("   - P1 snoops the bus and sees it has X in 'Exclusive' state.")
    print("   - P1 provides the data to P2 via a cache-to-cache transfer. (1 message)")
    print("   - P1's state changes from E to 'Shared' (S). P2 loads X in 'Shared' (S) state.")
    messages_this_step = 2
    step_messages.append(messages_this_step)
    message_count += messages_this_step
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print(f"   Messages exchanged in this step: {messages_this_step}")
    print(f"   Total messages so far: {message_count}")
    print_state()

    # --- Step 3: P1 writes X = 1 ---
    print("3. P1 writes X = 1")
    print("   - P1 has the block in 'Shared' (S) state. To write, it must gain exclusive ownership.")
    print("   - P1 sends an 'Invalidate' message on the bus to invalidate other copies. (1 message)")
    print("   - P2 receives the invalidate message and changes its state from S to 'Invalid' (I).")
    print("   - P1's state changes from S to 'Modified' (M) after writing.")
    messages_this_step = 1
    step_messages.append(messages_this_step)
    message_count += messages_this_step
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print(f"   Messages exchanged in this step: {messages_this_step}")
    print(f"   Total messages so far: {message_count}")
    print_state()

    # --- Step 4: P3 reads X ---
    print("4. P3 reads X")
    print("   - P3 has a cache miss.")
    print("   - P3 sends a 'BusRd' message to the bus. (1 message)")
    print("   - P1 snoops and finds it has the block in 'Modified' (M) state.")
    print("   - P1 intervenes, sending the data on the bus. This data is read by P3 and also updates main memory (write-back). (1 message)")
    print("   - P1's state changes from M to S. P3 loads X in 'Shared' (S) state.")
    messages_this_step = 2
    step_messages.append(messages_this_step)
    message_count += messages_this_step
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print(f"   Messages exchanged in this step: {messages_this_step}")
    print(f"   Total messages so far: {message_count}")
    print_state()

    # --- Step 5: P2 writes X = 2 ---
    print("5. P2 writes X = 2")
    print("   - P2 has a write miss (state is I).")
    print("   - P2 sends a 'BusRdX' (Read for Ownership) message. (1 message)")
    print("   - P1 and P3 snoop the message and invalidate their 'Shared' copies (S -> I).")
    print("   - Main memory (now up-to-date) provides the data to P2. (1 message)")
    print("   - P2's state becomes 'Modified' (M) after writing.")
    messages_this_step = 2
    step_messages.append(messages_this_step)
    message_count += messages_this_step
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print(f"   Messages exchanged in this step: {messages_this_step}")
    print(f"   Total messages so far: {message_count}")
    print_state()
    
    # --- Step 6: P4 reads X ---
    print("6. P4 reads X")
    print("   - P4 has a cache miss.")
    print("   - P4 sends a 'BusRd' message. (1 message)")
    print("   - P2 snoops and finds it has the block in 'Modified' (M) state.")
    print("   - P2 intervenes, sending data to P4 and memory (write-back). (1 message)")
    print("   - P2's state changes from M to S. P4 loads X in 'Shared' (S) state.")
    messages_this_step = 2
    step_messages.append(messages_this_step)
    message_count += messages_this_step
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print(f"   Messages exchanged in this step: {messages_this_step}")
    print(f"   Total messages so far: {message_count}")
    print_state()

    # --- Step 7: P1 reads X ---
    print("7. P1 reads X")
    print("   - P1 has a cache miss.")
    print("   - P1 sends a 'BusRd' message. (1 message)")
    print("   - P2 and P4 have the block in 'Shared' state. One of them (e.g., P2) responds with the data via cache-to-cache transfer. (1 message)")
    print("   - P1 loads the data in 'Shared' (S) state.")
    messages_this_step = 2
    step_messages.append(messages_this_step)
    message_count += messages_this_step
    caches['P1'] = 'S'
    print(f"   Messages exchanged in this step: {messages_this_step}")
    print(f"   Total messages so far: {message_count}")
    print_state()

    # --- Final Calculation ---
    print("\nFinal Calculation:")
    equation = " + ".join(map(str, step_messages))
    print(f"The total number of cache coherence messages is the sum of messages from each step:")
    print(f"{equation} = {message_count}")

if __name__ == '__main__':
    solve_mesi_coherence()