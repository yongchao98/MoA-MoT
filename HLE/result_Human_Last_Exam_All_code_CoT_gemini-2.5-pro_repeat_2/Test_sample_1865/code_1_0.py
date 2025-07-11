def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a given sequence of operations and
    counts the coherence messages.
    """
    # Initial state
    states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    total_messages = 0
    message_log = []

    print("Initial State:")
    print(f"Caches: {states}")
    print(f"Total Messages: {total_messages}\n")

    # --- Step 1: P1 reads X ---
    step = 1
    messages_this_step = 1
    print(f"--- Step {step}: P1 reads X ---")
    print("- P1 has a Read Miss, as its cache state for X is 'I'.")
    print("- P1 sends a 'Read Miss' message on the bus.")
    print("- No other cache has a copy, so main memory provides the data.")
    print("- P1's cache line for X transitions to 'Exclusive' (E).")
    states['P1'] = 'E'
    total_messages += messages_this_step
    message_log.append(str(messages_this_step))
    print(f"Messages this step: {messages_this_step}")
    print(f"Current States: {states}")
    print(f"Total Messages: {total_messages}\n")

    # --- Step 2: P2 reads X ---
    step = 2
    messages_this_step = 1
    print(f"--- Step {step}: P2 reads X ---")
    print("- P2 has a Read Miss (state is 'I').")
    print("- P2 sends a 'Read Miss' message on the bus.")
    print("- P1 snoops the bus and sees the request. Since P1 is in state 'E', it provides the data directly to P2.")
    print("- Both P1 and P2 transition to the 'Shared' (S) state.")
    states['P1'] = 'S'
    states['P2'] = 'S'
    total_messages += messages_this_step
    message_log.append(str(messages_this_step))
    print(f"Messages this step: {messages_this_step}")
    print(f"Current States: {states}")
    print(f"Total Messages: {total_messages}\n")

    # --- Step 3: P1 writes X = 1 ---
    step = 3
    messages_this_step = 1
    print(f"--- Step {step}: P1 writes X = 1 ---")
    print("- P1 needs to write, but its state is 'S'. It must gain exclusive ownership.")
    print("- P1 sends an 'Invalidate' message on the bus.")
    print("- P2 snoops the bus, receives the invalidate, and changes its state for X to 'Invalid' (I).")
    print("- P1 transitions to the 'Modified' (M) state and updates its value of X.")
    states['P1'] = 'M'
    states['P2'] = 'I'
    total_messages += messages_this_step
    message_log.append(str(messages_this_step))
    print(f"Messages this step: {messages_this_step}")
    print(f"Current States: {states}")
    print(f"Total Messages: {total_messages}\n")

    # --- Step 4: P3 reads X ---
    step = 4
    messages_this_step = 2
    print(f"--- Step {step}: P3 reads X ---")
    print("- P3 has a Read Miss (state is 'I').")
    print("- P3 sends a 'Read Miss' message. (Message 1)")
    print("- P1 snoops the bus and sees the request. Since P1 has the data in 'Modified' (M) state, it must intervene.")
    print("- P1 provides the data to P3 and also writes the data back to main memory. This Write-Back is a separate bus transaction. (Message 2)")
    print("- P1 transitions from 'M' to 'S'. P3 transitions to 'S'.")
    states['P1'] = 'S'
    states['P3'] = 'S'
    total_messages += messages_this_step
    message_log.append(str(messages_this_step))
    print(f"Messages this step: {messages_this_step} (Read Miss + Write-Back)")
    print(f"Current States: {states}")
    print(f"Total Messages: {total_messages}\n")

    # --- Step 5: P2 writes X = 2 ---
    step = 5
    messages_this_step = 1
    print(f"--- Step {step}: P2 writes X = 2 ---")
    print("- P2 has a Write Miss (state is 'I').")
    print("- P2 sends a 'Write Miss' (or Read-with-Intent-to-Modify) message on the bus.")
    print("- This message serves to both request the data and invalidate other copies.")
    print("- P1 and P3 snoop the bus and invalidate their copies ('S' -> 'I').")
    print("- P2 receives the data, updates it, and transitions to 'Modified' (M).")
    states['P1'] = 'I'
    states['P3'] = 'I'
    states['P2'] = 'M'
    total_messages += messages_this_step
    message_log.append(str(messages_this_step))
    print(f"Messages this step: {messages_this_step}")
    print(f"Current States: {states}")
    print(f"Total Messages: {total_messages}\n")

    # --- Step 6: P4 reads X ---
    step = 6
    messages_this_step = 2
    print(f"--- Step {step}: P4 reads X ---")
    print("- P4 has a Read Miss (state is 'I').")
    print("- P4 sends a 'Read Miss' message. (Message 1)")
    print("- P2 snoops the bus, and since it has the data in state 'M', it must respond.")
    print("- P2 provides the data to P4 and writes it back to memory. This Write-Back is the second message. (Message 2)")
    print("- P2 transitions from 'M' to 'S'. P4 transitions to 'S'.")
    states['P2'] = 'S'
    states['P4'] = 'S'
    total_messages += messages_this_step
    message_log.append(str(messages_this_step))
    print(f"Messages this step: {messages_this_step} (Read Miss + Write-Back)")
    print(f"Current States: {states}")
    print(f"Total Messages: {total_messages}\n")

    # --- Step 7: P1 reads X ---
    step = 7
    messages_this_step = 1
    print(f"--- Step {step}: P1 reads X ---")
    print("- P1 has a Read Miss (state is 'I').")
    print("- P1 sends a 'Read Miss' message on the bus.")
    print("- P2 and P4 have the data in 'Shared' (S) state. One of them (or memory) supplies the data.")
    print("- No write-back is needed. P1 transitions to 'S'.")
    states['P1'] = 'S'
    total_messages += messages_this_step
    message_log.append(str(messages_this_step))
    print(f"Messages this step: {messages_this_step}")
    print(f"Final States: {states}")
    print(f"Total Messages: {total_messages}\n")
    
    # --- Final Calculation ---
    print("--- Summary ---")
    final_equation = " + ".join(message_log)
    print(f"Total messages = {final_equation} = {total_messages}")

solve_mesi_coherence()
print("<<<9>>>")