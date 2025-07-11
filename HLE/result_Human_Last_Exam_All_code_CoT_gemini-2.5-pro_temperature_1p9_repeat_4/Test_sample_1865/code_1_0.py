def solve_mesi_messages():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    A "message" is a bus request like Read Miss, Invalidate, or Read-For-Ownership.
    """
    # Initial state
    # Caches: I=Invalid, E=Exclusive, S=Shared, M=Modified
    states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    memory_value = 0
    total_messages = 0
    messages_per_step = []

    print("--- MESI Cache Coherence Simulation ---\n")

    # Step-by-step simulation
    # Operation 1: P1 reads X
    step = 1
    messages_this_step = 1
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    states['P1'] = 'E'
    print(f"Step {step}: P1 reads X")
    print(f"  - P1 has a miss (I -> E), sends a Read Miss message.")
    print(f"  - Messages: {messages_this_step}. States: {states}\n")
    
    # Operation 2: P2 reads X
    step = 2
    messages_this_step = 1
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    states['P1'] = 'S'
    states['P2'] = 'S'
    print(f"Step {step}: P2 reads X")
    print(f"  - P2 has a miss, sends a Read Miss message.")
    print(f"  - P1 supplies data (E -> S). P2 enters S state.")
    print(f"  - Messages: {messages_this_step}. States: {states}\n")

    # Operation 3: P1 writes X = 1
    step = 3
    messages_this_step = 1
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    states['P1'] = 'M'
    states['P2'] = 'I'
    memory_value = 0 # Memory is now stale
    print(f"Step {step}: P1 writes X = 1")
    print(f"  - P1 must upgrade (S -> M), sends an Invalidate message.")
    print(f"  - P2 invalidates its copy (S -> I).")
    print(f"  - Messages: {messages_this_step}. States: {states}\n")
    
    # Operation 4: P3 reads X
    step = 4
    messages_this_step = 1
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    states['P1'] = 'S'
    states['P3'] = 'S'
    memory_value = 1 # P1 writes back to memory
    print(f"Step {step}: P3 reads X")
    print(f"  - P3 has a miss, sends a Read Miss message.")
    print(f"  - P1 supplies data (M -> S) and updates memory. P3 enters S state.")
    print(f"  - Messages: {messages_this_step}. States: {states}\n")

    # Operation 5: P2 writes X = 2
    step = 5
    messages_this_step = 1
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    states['P1'] = 'I'
    states['P2'] = 'M'
    states['P3'] = 'I'
    memory_value = 1 # Memory is now stale
    print(f"Step {step}: P2 writes X = 2")
    print(f"  - P2 has a write miss (I -> M), sends a Read-For-Ownership (RFO) message.")
    print(f"  - P1 and P3 invalidate their copies (S -> I).")
    print(f"  - Messages: {messages_this_step}. States: {states}\n")

    # Operation 6: P4 reads X
    step = 6
    messages_this_step = 1
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    states['P2'] = 'S'
    states['P4'] = 'S'
    memory_value = 2 # P2 writes back to memory
    print(f"Step {step}: P4 reads X")
    print(f"  - P4 has a miss, sends a Read Miss message.")
    print(f"  - P2 supplies data (M -> S) and updates memory. P4 enters S state.")
    print(f"  - Messages: {messages_this_step}. States: {states}\n")
    
    # Operation 7: P1 reads X
    step = 7
    messages_this_step = 1
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    states['P1'] = 'S'
    print(f"Step {step}: P1 reads X")
    print(f"  - P1 has a miss (I -> S), sends a Read Miss message.")
    print(f"  - A sharer (P2 or P4) supplies the data. P1 enters S state.")
    print(f"  - Messages: {messages_this_step}. States: {states}\n")

    # Final Calculation
    print("--- Final Calculation ---")
    equation = " + ".join(map(str, messages_per_step))
    print(f"Total messages = {equation} = {total_messages}")

solve_mesi_messages()