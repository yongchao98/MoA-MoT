def solve_mesi_problem():
    """
    Analyzes a sequence of processor operations under the MESI protocol
    and calculates the total number of cache coherence messages exchanged.
    """
    # Initialize the state of the system
    # Processor cache states for variable X: M(odified), E(xclusive), S(hared), I(nvalid)
    states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    
    # Value in memory
    memory_value = 0
    
    # Message counters
    total_messages = 0
    message_counts_per_step = []

    print("--- MESI Protocol Analysis ---")
    print(f"Initial State: P1({states['P1']}), P2({states['P2']}), P3({states['P3']}), P4({states['P4']})")
    print(f"Initial Memory Value of X: {memory_value}")
    print(f"Initial Total Messages: {total_messages}\n")

    # --- Step 1: P1 reads X ---
    messages_step = 2
    total_messages += messages_step
    message_counts_per_step.append(messages_step)
    print("1. P1 reads X")
    print("   - Action: P1 has a Read Miss (state is I).")
    print("   - Bus Activity:")
    print("     - P1 sends a 'Read Miss' request to the bus. (1 message)")
    print("     - Memory provides the data for X, as no cache has it. (1 message)")
    print("   - State Changes:")
    states['P1'] = 'E'  # Becomes Exclusive as it's the only cache with the data.
    print(f"     - P1 state: I -> E")
    print(f"   - Messages this step: {messages_step}. Total messages: {total_messages}\n")

    # --- Step 2: P2 reads X ---
    messages_step = 2
    total_messages += messages_step
    message_counts_per_step.append(messages_step)
    print("2. P2 reads X")
    print("   - Action: P2 has a Read Miss (state is I).")
    print("   - Bus Activity:")
    print("     - P2 sends a 'Read Miss' request to the bus. (1 message)")
    print("     - P1 snoops the request, has X in state E, and supplies the data (cache-to-cache transfer). (1 message)")
    print("   - State Changes:")
    states['P1'] = 'S'  # No longer exclusive.
    states['P2'] = 'S'  # Now shares the data.
    print(f"     - P1 state: E -> S")
    print(f"     - P2 state: I -> S")
    print(f"   - Messages this step: {messages_step}. Total messages: {total_messages}\n")

    # --- Step 3: P1 writes X = 1 ---
    messages_step = 1
    total_messages += messages_step
    message_counts_per_step.append(messages_step)
    print("3. P1 writes X = 1")
    print("   - Action: P1 has a Write Hit, but must upgrade to Modified state.")
    print("   - Bus Activity:")
    print("     - P1 sends an 'Invalidate' request on the bus to claim exclusive ownership. (1 message)")
    print("   - State Changes:")
    states['P1'] = 'M'  # Gains exclusive ownership and modifies data.
    states['P2'] = 'I'  # P2's copy is invalidated.
    print(f"     - P1 state: S -> M")
    print(f"     - P2 state: S -> I")
    print(f"   - Messages this step: {messages_step}. Total messages: {total_messages}\n")

    # --- Step 4: P3 reads X ---
    messages_step = 2
    total_messages += messages_step
    message_counts_per_step.append(messages_step)
    print("4. P3 reads X")
    print("   - Action: P3 has a Read Miss (state is I).")
    print("   - Bus Activity:")
    print("     - P3 sends a 'Read Miss' request to the bus. (1 message)")
    print("     - P1 snoops, finds it has X in M state, supplies the data, and writes back to memory. (1 message)")
    print("   - State Changes:")
    memory_value = 1  # P1 writes back the value 1.
    states['P1'] = 'S'  # Becomes shared after providing the data.
    states['P3'] = 'S'  # Receives the data.
    print(f"     - P1 state: M -> S")
    print(f"     - P3 state: I -> S")
    print("     - Main memory is updated to 1.")
    print(f"   - Messages this step: {messages_step}. Total messages: {total_messages}\n")

    # --- Step 5: P2 writes X = 2 ---
    messages_step = 2
    total_messages += messages_step
    message_counts_per_step.append(messages_step)
    print("5. P2 writes X = 2")
    print("   - Action: P2 has a Write Miss (state is I).")
    print("   - Bus Activity:")
    print("     - P2 sends a 'Read for Ownership' (RFO) request, which is a read and invalidate combined. (1 message)")
    print("     - P1 or P3 (in state S) supplies the data. (1 message)")
    print("   - State Changes:")
    states['P1'] = 'I'
    states['P3'] = 'I'
    states['P2'] = 'M'
    print(f"     - P1 and P3 states: S -> I")
    print(f"     - P2 state: I -> M")
    print(f"   - Messages this step: {messages_step}. Total messages: {total_messages}\n")

    # --- Step 6: P4 reads X ---
    messages_step = 2
    total_messages += messages_step
    message_counts_per_step.append(messages_step)
    print("6. P4 reads X")
    print("   - Action: P4 has a Read Miss (state is I).")
    print("   - Bus Activity:")
    print("     - P4 sends a 'Read Miss' request to the bus. (1 message)")
    print("     - P2 snoops, has X in M state, supplies the data, and writes back to memory. (1 message)")
    print("   - State Changes:")
    memory_value = 2  # P2 writes back the value 2.
    states['P2'] = 'S'
    states['P4'] = 'S'
    print(f"     - P2 state: M -> S")
    print(f"     - P4 state: I -> S")
    print("     - Main memory is updated to 2.")
    print(f"   - Messages this step: {messages_step}. Total messages: {total_messages}\n")

    # --- Step 7: P1 reads X ---
    messages_step = 2
    total_messages += messages_step
    message_counts_per_step.append(messages_step)
    print("7. P1 reads X")
    print("   - Action: P1 has a Read Miss (state is I).")
    print("   - Bus Activity:")
    print("     - P1 sends a 'Read Miss' request to the bus. (1 message)")
    print("     - P2 or P4 (in state S) supplies the data via cache-to-cache transfer. (1 message)")
    print("   - State Changes:")
    states['P1'] = 'S'
    print(f"     - P1 state: I -> S")
    print(f"   - Messages this step: {messages_step}. Total messages: {total_messages}\n")
    
    # --- Final Calculation ---
    print("--- Summary of Messages ---")
    equation = " + ".join(map(str, message_counts_per_step))
    print(f"Final Calculation: {equation} = {total_messages}")

if __name__ == '__main__':
    solve_mesi_problem()