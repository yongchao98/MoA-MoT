def solve_mesi_messages():
    """
    Simulates the MESI protocol for a sequence of operations and counts the messages.
    """
    # Initial states
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    total_messages = 0
    messages_per_step = []

    print("--- MESI Cache Coherence Simulation ---\n")
    print(f"Initial State: Caches = {caches}, Total Messages = {total_messages}\n" + "="*45)

    # 1. P1 reads X
    step = 1
    messages_this_step = 2  # Read Miss + Data Response from Memory
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    caches['P1'] = 'E'
    print(f"Step {step}: P1 reads X")
    print(f"  - P1 sends Read Miss. Memory responds. P1 state becomes Exclusive (E).")
    print(f"  - Messages: {messages_this_step}. Total Messages: {total_messages}")
    print(f"  - Caches: {caches}\n" + "="*45)

    # 2. P2 reads X
    step = 2
    messages_this_step = 2  # Read Miss + Data Response from P1's cache
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print(f"Step {step}: P2 reads X")
    print(f"  - P2 sends Read Miss. P1 provides data. P1 and P2 states become Shared (S).")
    print(f"  - Messages: {messages_this_step}. Total Messages: {total_messages}")
    print(f"  - Caches: {caches}\n" + "="*45)

    # 3. P1 writes X = 1
    step = 3
    messages_this_step = 1  # Invalidate message
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print(f"Step {step}: P1 writes X = 1")
    print(f"  - P1 sends Invalidate. P2 invalidates its copy. P1 state becomes Modified (M).")
    print(f"  - Messages: {messages_this_step}. Total Messages: {total_messages}")
    print(f"  - Caches: {caches}\n" + "="*45)

    # 4. P3 reads X
    step = 4
    messages_this_step = 2  # Read Miss + Data Response from P1
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print(f"Step {step}: P3 reads X")
    print(f"  - P3 sends Read Miss. P1 provides data. P1 and P3 states become Shared (S).")
    print(f"  - Messages: {messages_this_step}. Total Messages: {total_messages}")
    print(f"  - Caches: {caches}\n" + "="*45)

    # 5. P2 writes X = 2
    step = 5
    messages_this_step = 2  # Read For Ownership (RFO) + Data Response
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    caches['P1'] = 'I'
    caches['P2'] = 'M'
    caches['P3'] = 'I'
    print(f"Step {step}: P2 writes X = 2")
    print(f"  - P2 sends RFO. P1 provides data. P1/P3 invalidate. P2 state becomes Modified (M).")
    print(f"  - Messages: {messages_this_step}. Total Messages: {total_messages}")
    print(f"  - Caches: {caches}\n" + "="*45)

    # 6. P4 reads X
    step = 6
    messages_this_step = 2  # Read Miss + Data Response from P2
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print(f"Step {step}: P4 reads X")
    print(f"  - P4 sends Read Miss. P2 provides data. P2 and P4 states become Shared (S).")
    print(f"  - Messages: {messages_this_step}. Total Messages: {total_messages}")
    print(f"  - Caches: {caches}\n" + "="*45)

    # 7. P1 reads X
    step = 7
    messages_this_step = 2  # Read Miss + Data Response from P2/P4
    total_messages += messages_this_step
    messages_per_step.append(messages_this_step)
    caches['P1'] = 'S'
    print(f"Step {step}: P1 reads X")
    print(f"  - P1 sends Read Miss. P2/P4 provides data. P1 state becomes Shared (S).")
    print(f"  - Messages: {messages_this_step}. Total Messages: {total_messages}")
    print(f"  - Caches: {caches}\n" + "="*45)

    # Final Result
    equation = " + ".join(map(str, messages_per_step))
    print("\n--- Final Calculation ---")
    print(f"Total messages exchanged is the sum of messages from each step:")
    print(f"{equation} = {total_messages}")

solve_mesi_messages()
<<<13>>>