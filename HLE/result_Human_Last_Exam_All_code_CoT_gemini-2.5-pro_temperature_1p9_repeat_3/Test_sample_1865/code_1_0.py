def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a given sequence of operations
    and calculates the number of cache coherence messages.
    """
    # Caches stores the state of variable X for each processor
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    message_count = 0
    
    # Store messages per step for the final equation
    step_messages = []

    print("Initial State:")
    print(f"Caches: {caches}")
    print(f"Total Messages: {message_count}\n" + "="*40)

    # --- Step 1: P1 reads X ---
    # P1 has a read miss, issues a BusRd (Read Miss) message.
    # Memory responds. P1's state becomes Exclusive (E) as it's the only cache with X.
    messages_in_step = 1
    message_count += messages_in_step
    step_messages.append(messages_in_step)
    caches['P1'] = 'E'
    print("Step 1: P1 reads X")
    print("Action: P1 issues a Read Miss (BusRd) message.")
    print("State Change: P1 -> Exclusive (E).")
    print(f"Caches: {caches}")
    print(f"Messages for this step: {messages_in_step}")
    print(f"Total Messages: {message_count}\n" + "="*40)

    # --- Step 2: P2 reads X ---
    # P2 has a read miss, issues a BusRd.
    # P1 snoops, has X in E state, provides the data, and changes its state to Shared (S).
    # P2 loads data and sets its state to S.
    messages_in_step = 1
    message_count += messages_in_step
    step_messages.append(messages_in_step)
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print("Step 2: P2 reads X")
    print("Action: P2 issues a Read Miss (BusRd) message. P1 supplies the data.")
    print("State Change: P1 -> Shared (S), P2 -> Shared (S).")
    print(f"Caches: {caches}")
    print(f"Messages for this step: {messages_in_step}")
    print(f"Total Messages: {message_count}\n" + "="*40)

    # --- Step 3: P1 writes X = 1 ---
    # P1 wants to write to a line in S state. It must invalidate other copies.
    # P1 issues a BusUpgr (Invalidate) message.
    # P2 invalidates its copy (S -> I). P1's state becomes Modified (M).
    messages_in_step = 1
    message_count += messages_in_step
    step_messages.append(messages_in_step)
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print("Step 3: P1 writes X = 1")
    print("Action: P1 issues an Invalidate (BusUpgr) message.")
    print("State Change: P1 -> Modified (M), P2 -> Invalid (I).")
    print(f"Caches: {caches}")
    print(f"Messages for this step: {messages_in_step}")
    print(f"Total Messages: {message_count}\n" + "="*40)

    # --- Step 4: P3 reads X ---
    # P3 has a read miss, issues a BusRd.
    # P1 snoops, has X in M (dirty) state. It must write back to memory AND provide data to P3.
    # This costs two messages: P3's BusRd and P1's Write-Back.
    # P1's state becomes S. P3's state becomes S.
    messages_in_step = 2
    message_count += messages_in_step
    step_messages.append(messages_in_step)
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print("Step 4: P3 reads X")
    print("Action: P3 issues a Read Miss (BusRd). P1 has a dirty copy, so it writes back to memory and supplies data.")
    print("State Change: P1 -> Shared (S), P3 -> Shared (S).")
    print(f"Caches: {caches}")
    print(f"Messages for this step: {messages_in_step} (1 Read Miss + 1 Write-Back)")
    print(f"Total Messages: {message_count}\n" + "="*40)
    
    # --- Step 5: P2 writes X = 2 ---
    # P2 has a write miss (state is I), issues a BusRdX (Read-For-Ownership) message.
    # This message also invalidates other copies. P1 and P3 go from S to I.
    # P2 loads data and its state becomes M.
    messages_in_step = 1
    message_count += messages_in_step
    step_messages.append(messages_in_step)
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print("Step 5: P2 writes X = 2")
    print("Action: P2 issues a Read-For-Ownership (BusRdX) message.")
    print("State Change: P1 -> Invalid (I), P2 -> Modified (M), P3 -> Invalid (I).")
    print(f"Caches: {caches}")
    print(f"Messages for this step: {messages_in_step}")
    print(f"Total Messages: {message_count}\n" + "="*40)

    # --- Step 6: P4 reads X ---
    # P4 has a read miss, issues a BusRd.
    # P2 snoops, has X in M state. It performs a Write-Back and provides data to P4.
    # This costs two messages: P4's BusRd and P2's Write-Back.
    # P2's state becomes S. P4's state becomes S.
    messages_in_step = 2
    message_count += messages_in_step
    step_messages.append(messages_in_step)
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print("Step 6: P4 reads X")
    print("Action: P4 issues a Read Miss (BusRd). P2 has a dirty copy, so it writes back to memory and supplies data.")
    print("State Change: P2 -> Shared (S), P4 -> Shared (S).")
    print(f"Caches: {caches}")
    print(f"Messages for this step: {messages_in_step} (1 Read Miss + 1 Write-Back)")
    print(f"Total Messages: {message_count}\n" + "="*40)

    # --- Step 7: P1 reads X ---
    # P1 has a read miss, issues a BusRd.
    # Data is available in S state from P2/P4 (and also memory). One of them supplies the data.
    # P1's state becomes S. No change for P2/P4.
    messages_in_step = 1
    message_count += messages_in_step
    step_messages.append(messages_in_step)
    caches['P1'] = 'S'
    print("Step 7: P1 reads X")
    print("Action: P1 issues a Read Miss (BusRd). Data is supplied by another cache or memory.")
    print("State Change: P1 -> Shared (S).")
    print(f"Caches: {caches}")
    print(f"Messages for this step: {messages_in_step}")
    print(f"Total Messages: {message_count}\n" + "="*40)

    # --- Final Calculation ---
    equation = " + ".join(map(str, step_messages))
    print("Final Calculation:")
    print(f"Total Messages = {equation} = {sum(step_messages)}")

solve_mesi_coherence()
<<<9>>>