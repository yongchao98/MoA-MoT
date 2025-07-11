def solve_mesi_messages():
    """
    Traces the MESI protocol for a given sequence of operations and
    calculates the total number of cache coherence messages.
    """
    # M: Modified, E: Exclusive, S: Shared, I: Invalid
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    total_messages = 0
    messages_per_step = []

    print("Initial State:")
    print(f"  Caches: {caches}")
    print(f"  Total Messages: {total_messages}\n")

    # Step 1: P1 reads X
    step = 1
    step_messages = 1
    messages_per_step.append(step_messages)
    total_messages += step_messages
    caches['P1'] = 'E'
    print(f"Step {step}: P1 reads X")
    print(f"  - P1 has a read miss, issues a 'BusRd' message ({step_messages} message).")
    print("  - Memory responds. P1 is the only cache with X.")
    print(f"  - Cache states change: P1 -> E")
    print(f"  - Current Caches: {caches}")
    print(f"  - Running Total Messages: {total_messages}\n")

    # Step 2: P2 reads X
    step = 2
    step_messages = 1
    messages_per_step.append(step_messages)
    total_messages += step_messages
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print(f"Step {step}: P2 reads X")
    print(f"  - P2 has a read miss, issues a 'BusRd' message ({step_messages} message).")
    print("  - P1 snoops, supplies data to P2, and changes its own state to S.")
    print(f"  - Cache states change: P1 -> S, P2 -> S")
    print(f"  - Current Caches: {caches}")
    print(f"  - Running Total Messages: {total_messages}\n")

    # Step 3: P1 writes X = 1
    step = 3
    step_messages = 1
    messages_per_step.append(step_messages)
    total_messages += step_messages
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print(f"Step {step}: P1 writes X = 1")
    print(f"  - P1 needs to invalidate other copies to write. Issues a 'BusRdX' (Invalidate) message ({step_messages} message).")
    print("  - P2 snoops and invalidates its copy.")
    print(f"  - Cache states change: P1 -> M, P2 -> I")
    print(f"  - Current Caches: {caches}")
    print(f"  - Running Total Messages: {total_messages}\n")

    # Step 4: P3 reads X
    step = 4
    step_messages = 2
    messages_per_step.append(step_messages)
    total_messages += step_messages
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print(f"Step {step}: P3 reads X")
    print("  - P3 has a read miss, issues a 'BusRd' message (1 message).")
    print("  - P1 snoops, sees it has a Modified copy, so it must first send a 'Write-Back' message to update memory (1 message).")
    print(f"  - Total messages for step: {step_messages}.")
    print(f"  - Cache states change: P1 -> S, P3 -> S")
    print(f"  - Current Caches: {caches}")
    print(f"  - Running Total Messages: {total_messages}\n")

    # Step 5: P2 writes X = 2
    step = 5
    step_messages = 1
    messages_per_step.append(step_messages)
    total_messages += step_messages
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print(f"Step {step}: P2 writes X = 2")
    print(f"  - P2 has a write miss, issues a 'BusRdX' message to get the data and invalidate other copies ({step_messages} message).")
    print("  - P1 and P3 snoop and invalidate their copies.")
    print(f"  - Cache states change: P2 -> M, P1 -> I, P3 -> I")
    print(f"  - Current Caches: {caches}")
    print(f"  - Running Total Messages: {total_messages}\n")

    # Step 6: P4 reads X
    step = 6
    step_messages = 2
    messages_per_step.append(step_messages)
    total_messages += step_messages
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print(f"Step {step}: P4 reads X")
    print("  - P4 has a read miss, issues a 'BusRd' message (1 message).")
    print("  - P2 snoops, has a Modified copy, sends a 'Write-Back' message to memory (1 message).")
    print(f"  - Total messages for step: {step_messages}.")
    print(f"  - Cache states change: P2 -> S, P4 -> S")
    print(f"  - Current Caches: {caches}")
    print(f"  - Running Total Messages: {total_messages}\n")
    
    # Step 7: P1 reads X
    step = 7
    step_messages = 1
    messages_per_step.append(step_messages)
    total_messages += step_messages
    caches['P1'] = 'S'
    print(f"Step {step}: P1 reads X")
    print(f"  - P1 has a read miss, issues a 'BusRd' message ({step_messages} message).")
    print("  - P2 (or P4) has a Shared copy and provides the data. No Write-Back is needed.")
    print(f"  - Cache states change: P1 -> S")
    print(f"  - Current Caches: {caches}")
    print(f"  - Running Total Messages: {total_messages}\n")
    
    # Final Calculation
    print("---")
    print("Final Calculation:")
    equation_str = " + ".join(map(str, messages_per_step))
    print(f"Total Messages = {equation_str} = {total_messages}")

solve_mesi_messages()
<<<9>>>