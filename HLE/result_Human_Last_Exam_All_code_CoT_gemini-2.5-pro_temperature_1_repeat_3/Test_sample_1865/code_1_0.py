def solve_mesi_messages():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # Processor caches, initialized to Invalid
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    
    # Sequence of operations: (processor, 'read' or 'write')
    operations = [
        ('P1', 'read'),
        ('P2', 'read'),
        ('P1', 'write'),
        ('P3', 'read'),
        ('P2', 'write'),
        ('P4', 'read'),
        ('P1', 'read'),
    ]

    operation_descriptions = [
        "P1 reads X", "P2 reads X", "P1 writes X = 1", "P3 reads X",
        "P2 writes X = 2", "P4 reads X", "P1 reads X"
    ]

    total_messages = 0
    message_counts_per_step = []

    print(f"Initial State: {caches}\n")

    for i, (proc, op_type) in enumerate(operations):
        messages_this_step = 0
        print(f"Step {i+1}: {operation_descriptions[i]}")
        
        current_state = caches[proc]
        other_procs = [p for p in caches if p != proc]
        sharers = [p for p, s in caches.items() if s == 'S']
        exclusive_holder = next((p for p, s in caches.items() if s in ['E', 'M']), None)

        if op_type == 'read':
            if current_state == 'I':
                # Read Miss
                messages_this_step = 1
                print("  - Action: Read Miss. Issues a bus message (Read Miss).")
                if exclusive_holder: # Another cache has it in E or M
                    if caches[exclusive_holder] == 'E':
                         print(f"  - Response: {exclusive_holder} (in E) supplies data. Both {proc} and {exclusive_holder} become S.")
                         caches[exclusive_holder] = 'S'
                    else: # State is M
                         print(f"  - Response: {exclusive_holder} (in M) supplies data and writes to memory. Both {proc} and {exclusive_holder} become S.")
                         caches[exclusive_holder] = 'S'
                    caches[proc] = 'S'
                elif sharers: # Other caches have it in S
                    print(f"  - Response: A sharer (e.g., {sharers[0]}) supplies data. {proc} becomes S.")
                    caches[proc] = 'S'
                else: # No other cache has the data
                    print(f"  - Response: Memory supplies data. {proc} becomes E.")
                    caches[proc] = 'E'
        
        elif op_type == 'write':
            if current_state in ['I', 'S']:
                # Write Miss (from I) or Write Hit on Shared (from S)
                # Both require a bus message to invalidate other copies.
                messages_this_step = 1
                if current_state == 'I':
                    print("  - Action: Write Miss. Issues a bus message (Read-Exclusive).")
                else: # current_state == 'S'
                    print("  - Action: Write Hit on Shared. Issues a bus message (Write Invalidate).")

                # Invalidate all other copies
                for p in other_procs:
                    if caches[p] != 'I':
                        print(f"  - Response: {p} invalidates its copy (S/E -> I).")
                        caches[p] = 'I'
                
                caches[proc] = 'M'
                print(f"  - Result: {proc} becomes M.")
        
        total_messages += messages_this_step
        message_counts_per_step.append(str(messages_this_step))
        print(f"  - Current Messages: {total_messages}")
        print(f"  - New States: {caches}\n")

    # Format the final equation string
    equation = " + ".join(message_counts_per_step)
    print(f"Total messages = {equation} = {total_messages}")


solve_mesi_messages()
<<<7>>>