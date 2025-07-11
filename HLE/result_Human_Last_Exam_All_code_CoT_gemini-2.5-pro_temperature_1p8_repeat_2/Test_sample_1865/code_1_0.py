def solve_mesi_messages():
    """
    This function simulates the MESI protocol for a given sequence of operations
    and calculates the total number of cache coherence messages.
    """
    # Initial states for X in each processor's cache: I = Invalid
    cache_states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    total_messages = 0
    messages_per_step = []

    print("Analyzing the sequence of operations with the MESI protocol...\n")

    # 1. P1 reads X
    operation = "1. P1 reads X"
    messages_this_step = 1 # P1 Read Miss (BusRd)
    cache_states['P1'] = 'E'
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    print(f"{operation}: P1 sends a 'Read Miss'. States: {cache_states}. Messages so far: {total_messages}")

    # 2. P2 reads X
    operation = "2. P2 reads X"
    messages_this_step = 1 # P2 Read Miss (BusRd)
    cache_states['P1'] = 'S'
    cache_states['P2'] = 'S'
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    print(f"{operation}: P2 sends a 'Read Miss'. States: {cache_states}. Messages so far: {total_messages}")

    # 3. P1 writes X = 1
    operation = "3. P1 writes X = 1"
    messages_this_step = 1 # P1 Invalidate (BusUpgr)
    cache_states['P1'] = 'M'
    cache_states['P2'] = 'I'
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    print(f"{operation}: P1 sends an 'Invalidate'. States: {cache_states}. Messages so far: {total_messages}")

    # 4. P3 reads X
    operation = "4. P3 reads X"
    messages_this_step = 1 # P3 Read Miss (BusRd)
    cache_states['P1'] = 'S'
    cache_states['P3'] = 'S'
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    print(f"{operation}: P3 sends a 'Read Miss'. States: {cache_states}. Messages so far: {total_messages}")

    # 5. P2 writes X = 2
    operation = "5. P2 writes X = 2"
    messages_this_step = 1 # P2 Write Miss -> Read Invalidate (BusRdX)
    cache_states['P1'] = 'I'
    cache_states['P2'] = 'M'
    cache_states['P3'] = 'I'
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    print(f"{operation}: P2 sends a 'Read Invalidate'. States: {cache_states}. Messages so far: {total_messages}")

    # 6. P4 reads X
    operation = "6. P4 reads X"
    messages_this_step = 1 # P4 Read Miss (BusRd)
    cache_states['P2'] = 'S'
    cache_states['P4'] = 'S'
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    print(f"{operation}: P4 sends a 'Read Miss'. States: {cache_states}. Messages so far: {total_messages}")

    # 7. P1 reads X
    operation = "7. P1 reads X"
    messages_this_step = 1 # P1 Read Miss (BusRd)
    cache_states['P1'] = 'S'
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    print(f"{operation}: P1 sends a 'Read Miss'. States: {cache_states}. Messages so far: {total_messages}")

    # Final Calculation
    equation = " + ".join(map(str, messages_per_step))
    print("\n--- Final Calculation ---")
    print(f"Total number of cache coherence messages = {equation} = {total_messages}")

solve_mesi_messages()
<<<7>>>