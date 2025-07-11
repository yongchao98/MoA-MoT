def solve_mesi_coherence():
    """
    Traces the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # Initial state
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    total_messages = 0
    step_messages = []

    def print_status(operation, messages_this_step, explanation):
        nonlocal total_messages
        total_messages += messages_this_step
        step_messages.append(messages_this_step)
        print(f"--- {operation} ---")
        print(f"Action: {explanation}")
        print(f"Messages: {messages_this_step}")
        print(f"New States: {caches}")
        print(f"Total Messages: {total_messages}\n")

    # --- Initial State ---
    print("Initial State: All caches are Invalid (I).\n")

    # 1. P1 reads X
    # P1 has a Read Miss. It issues a Read Miss on the bus.
    # No other cache has X, so memory responds with the data. P1's cache state becomes Exclusive (E).
    # Messages: 1 (Read Miss) + 1 (Data Response) = 2
    caches['P1'] = 'E'
    print_status("1. P1 reads X", 2, "Read Miss. Memory responds. P1->E.")
    
    # 2. P2 reads X
    # P2 has a Read Miss. It issues a Read Miss on the bus.
    # P1 snoops the bus and sees it has the line in state E. It supplies the data to P2.
    # Both P1 and P2 now have the data, so their states become Shared (S).
    # Messages: 1 (Read Miss) + 1 (Cache-to-Cache Data) = 2
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print_status("2. P2 reads X", 2, "Read Miss. P1 supplies data. P1->S, P2->S.")

    # 3. P1 writes X = 1
    # P1 has the line in state S. To write, it must invalidate other shared copies.
    # P1 issues an Invalidate message on the bus. P2 receives it, invalidates its copy (S->I), and sends an Acknowledgment.
    # P1's state becomes Modified (M).
    # Messages: 1 (Invalidate) + 1 (Invalidate Ack from P2) = 2
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print_status("3. P1 writes X = 1", 2, "Write Hit on Shared line. P1 sends Invalidate. P2 acknowledges. P1->M, P2->I.")

    # 4. P3 reads X
    # P3 has a Read Miss. It issues a Read Miss on the bus.
    # P1 snoops and sees it has the line in state M (dirty). It must supply the data to P3 and write it back to memory.
    # Both P1 and P3 states become Shared (S).
    # Messages: 1 (Read Miss) + 1 (Data supply from P1 to P3/memory) = 2
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print_status("4. P3 reads X", 2, "Read Miss. P1 has modified data, supplies it, and writes back. P1->S, P3->S.")

    # 5. P2 writes X = 2
    # P2 has the line in state I (a Write Miss).
    # P2 issues a Read-Invalidate (RWITM) message. P1 and P3 invalidate their copies (S->I) and send Acks.
    # Memory (now up-to-date) supplies the data to P2. P2's state becomes Modified (M).
    # Messages: 1 (Read-Invalidate) + 2 (Acks from P1, P3) + 1 (Data Response) = 4
    caches['P2'] = 'M'
    caches['P1'] = 'I'
    caches['P3'] = 'I'
    print_status("5. P2 writes X = 2", 4, "Write Miss. P2 sends Read-Invalidate. P1,P3 ack and invalidate. P2->M.")

    # 6. P4 reads X
    # P4 has a Read Miss. It issues a Read Miss on the bus.
    # P2 snoops, sees it has the line in state M. It supplies the data to P4 and writes it back to memory.
    # Both P2 and P4 states become Shared (S).
    # Messages: 1 (Read Miss) + 1 (Data supply from P2 to P4/memory) = 2
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print_status("6. P4 reads X", 2, "Read Miss. P2 has modified data, supplies it, and writes back. P2->S, P4->S.")

    # 7. P1 reads X
    # P1 has a Read Miss. It issues a Read Miss on the bus.
    # P2 and P4 have the data in state S. One of them (e.g., P2) supplies the data via a cache-to-cache transfer.
    # P1's state becomes Shared (S).
    # Messages: 1 (Read Miss) + 1 (Cache-to-Cache Data) = 2
    caches['P1'] = 'S'
    print_status("7. P1 reads X", 2, "Read Miss. P2 (or P4) supplies data. P1->S.")

    # --- Final Calculation ---
    print("--- Final Calculation ---")
    equation = " + ".join(map(str, step_messages))
    print(f"The total number of messages is the sum from each step:")
    print(f"{equation} = {total_messages}")

solve_mesi_coherence()
print("<<<16>>>")