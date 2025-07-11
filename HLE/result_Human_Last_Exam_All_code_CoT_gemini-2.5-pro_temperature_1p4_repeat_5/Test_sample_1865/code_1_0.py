def solve_mesi_coherence():
    """
    Calculates the number of cache coherence messages for a given sequence of operations
    in a 4-processor system using the MESI protocol.
    """
    
    # M: Modified, E: Exclusive, S: Shared, I: Invalid
    states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    memory_value = 0
    total_messages = 0
    messages_per_step = []

    print("Initial State: Variable X is not cached. All caches are in 'Invalid' state.\n")

    # Helper function to print current states
    def print_states(step, operation):
        print(f"Step {step}: {operation}")
        print(f"  - Cache States: P1={states['P1']}, P2={states['P2']}, P3={states['P3']}, P4={states['P4']}")

    # --- Step 1: P1 reads X ---
    step = 1
    operation = "P1 reads X"
    print_states(step, operation)
    # P1 has a read miss. It sends a BusRd message. Memory provides the data.
    # Since no other cache has the data, P1's cache state becomes Exclusive (E).
    messages = 1
    total_messages += messages
    messages_per_step.append(messages)
    states['P1'] = 'E'
    print(f"  - Action: P1 read miss. P1 issues a BusRd message.")
    print(f"  - Messages: {messages} (BusRd)")
    print(f"  - New States: P1->E\n")

    # --- Step 2: P2 reads X ---
    step = 2
    operation = "P2 reads X"
    print_states(step, operation)
    # P2 has a read miss. It sends a BusRd message.
    # P1 snoops the bus, sees the request, and provides the data (cache-to-cache transfer).
    # P1's state changes from Exclusive to Shared (S). P2's state becomes Shared (S).
    messages = 1
    total_messages += messages
    messages_per_step.append(messages)
    states['P1'] = 'S'
    states['P2'] = 'S'
    print(f"  - Action: P2 read miss. P2 issues a BusRd message. P1 supplies data.")
    print(f"  - Messages: {messages} (BusRd)")
    print(f"  - New States: P1->S, P2->S\n")

    # --- Step 3: P1 writes X = 1 ---
    step = 3
    operation = "P1 writes X = 1"
    print_states(step, operation)
    # P1 has the data in Shared state. To write, it must upgrade its state to Modified (M).
    # P1 sends an Invalidate message on the bus. P2 receives it and invalidates its copy.
    messages = 1
    total_messages += messages
    messages_per_step.append(messages)
    states['P1'] = 'M'
    states['P2'] = 'I'
    memory_value = 1
    print(f"  - Action: P1 write hit (on Shared data). P1 issues an Invalidate message.")
    print(f"  - Messages: {messages} (Invalidate)")
    print(f"  - New States: P1->M, P2->I\n")
    
    # --- Step 4: P3 reads X ---
    step = 4
    operation = "P3 reads X"
    print_states(step, operation)
    # P3 has a read miss and sends a BusRd message.
    # P1 has the data in Modified state. P1 provides the data to P3 and also writes it back to memory.
    # This involves two bus operations: P3's request and P1's write-back.
    # Both P1 and P3 end up in the Shared state.
    messages = 2
    total_messages += messages
    messages_per_step.append(messages)
    states['P1'] = 'S'
    states['P3'] = 'S'
    print(f"  - Action: P3 read miss. P3 issues BusRd. P1 has modified data, supplies it, and writes back to memory.")
    print(f"  - Messages: {messages} (BusRd + Write-Back)")
    print(f"  - New States: P1->S, P3->S\n")

    # --- Step 5: P2 writes X = 2 ---
    step = 5
    operation = "P2 writes X = 2"
    print_states(step, operation)
    # P2's copy is Invalid (a write miss). It needs exclusive ownership.
    # P2 sends a BusRdX (Read for Ownership) message.
    # P1 and P3 snoop this message and invalidate their copies. P2 gets the data, writes it, and sets its state to Modified.
    messages = 1
    total_messages += messages
    messages_per_step.append(messages)
    states['P1'] = 'I'
    states['P2'] = 'M'
    states['P3'] = 'I'
    memory_value = 2
    print(f"  - Action: P2 write miss. P2 issues a BusRdX message to get data and invalidate others.")
    print(f"  - Messages: {messages} (BusRdX)")
    print(f"  - New States: P1->I, P2->M, P3->I\n")

    # --- Step 6: P4 reads X ---
    step = 6
    operation = "P4 reads X"
    print_states(step, operation)
    # P4 has a read miss and sends a BusRd message.
    # P2 has the data in Modified state. It provides the data to P4 and writes it back to memory.
    # Similar to step 4, this is counted as two messages.
    messages = 2
    total_messages += messages
    messages_per_step.append(messages)
    states['P2'] = 'S'
    states['P4'] = 'S'
    print(f"  - Action: P4 read miss. P4 issues BusRd. P2 has modified data, supplies it, and writes back to memory.")
    print(f"  - Messages: {messages} (BusRd + Write-Back)")
    print(f"  - New States: P2->S, P4->S\n")

    # --- Step 7: P1 reads X ---
    step = 7
    operation = "P1 reads X"
    print_states(step, operation)
    # P1's copy is Invalid (a read miss). It sends a BusRd message.
    # P2 and P4 have the data in Shared state. One of them provides the data to P1.
    # All three (P1, P2, P4) will now have the data in Shared state.
    messages = 1
    total_messages += messages
    messages_per_step.append(messages)
    states['P1'] = 'S'
    print(f"  - Action: P1 read miss. P1 issues a BusRd message. P2 or P4 supplies data.")
    print(f"  - Messages: {messages} (BusRd)")
    print(f"  - New States: P1->S\n")

    # --- Final Calculation ---
    print("-----------------------------------------")
    print("Total Cache Coherence Messages Calculation:")
    equation = " + ".join(map(str, messages_per_step))
    print(f"Messages = {equation} = {total_messages}")
    print("-----------------------------------------")
    
    return total_messages

if __name__ == '__main__':
    final_answer = solve_mesi_coherence()
    print(f"<<<{final_answer}>>>")