def solve_mesi_messages():
    """
    Simulates the MESI protocol for a sequence of operations and counts the messages.
    """
    # Each processor's cache state for variable X. States: M, E, S, I
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    
    # List to store the number of messages for each step
    message_counts = []

    # --- Operation 1: P1 reads X ---
    # P1 read miss. Issues BusRd. No other cache has it. Fetches from memory.
    messages_step_1 = 1
    caches['P1'] = 'E'
    message_counts.append(messages_step_1)

    # --- Operation 2: P2 reads X ---
    # P2 read miss. Issues BusRd. P1 snoops, provides data.
    messages_step_2 = 1
    caches['P1'] = 'S' # P1 downgrades to Shared
    caches['P2'] = 'S' # P2 loads as Shared
    message_counts.append(messages_step_2)

    # --- Operation 3: P1 writes X = 1 ---
    # P1 has X in Shared state. Issues Invalidate to upgrade.
    messages_step_3 = 1
    caches['P1'] = 'M' # P1 becomes Modified
    if caches['P2'] == 'S': # P2 snoops and invalidates its copy
        caches['P2'] = 'I'
    message_counts.append(messages_step_3)

    # --- Operation 4: P3 reads X ---
    # P3 read miss. Issues BusRd. P1 (in Modified) snoops and provides data.
    messages_step_4 = 1
    caches['P1'] = 'S' # P1 downgrades to Shared
    caches['P3'] = 'S' # P3 loads as Shared
    message_counts.append(messages_step_4)
    
    # --- Operation 5: P2 writes X = 2 ---
    # P2 write miss (state is Invalid). Issues BusRdX.
    messages_step_5 = 1
    caches['P2'] = 'M' # P2 becomes Modified
    # All other caches snoop and invalidate their copies.
    if caches['P1'] == 'S':
        caches['P1'] = 'I'
    if caches['P3'] == 'S':
        caches['P3'] = 'I'
    message_counts.append(messages_step_5)

    # --- Operation 6: P4 reads X ---
    # P4 read miss. Issues BusRd. P2 (in Modified) snoops and provides data.
    messages_step_6 = 1
    caches['P2'] = 'S' # P2 downgrades to Shared
    caches['P4'] = 'S' # P4 loads as Shared
    message_counts.append(messages_step_6)
    
    # --- Operation 7: P1 reads X ---
    # P1 read miss (state is Invalid). Issues BusRd.
    messages_step_7 = 1
    # P2 or P4 provides data. P1 becomes Shared. Others remain Shared.
    caches['P1'] = 'S'
    message_counts.append(messages_step_7)

    # Calculate and print the final result
    total_messages = sum(message_counts)
    equation = " + ".join(map(str, message_counts))
    print(f"The number of cache coherence messages for each step are: {equation}")
    print(f"Total Messages = {equation} = {total_messages}")

solve_mesi_messages()
<<<7>>>