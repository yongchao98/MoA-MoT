def simulate_mesi_protocol():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # Cache states: I (Invalid), S (Shared), E (Exclusive), M (Modified)
    caches = {
        'P1': {'state': 'I', 'value': None},
        'P2': {'state': 'I', 'value': None},
        'P3': {'state': 'I', 'value': None},
        'P4': {'state': 'I', 'value': None},
    }
    memory_value = 0
    message_count = 0
    
    print("Initial State:")
    print(f"Memory Value: {memory_value}")
    print(f"Cache States: { {p: c['state'] for p, c in caches.items()} }")
    print(f"Total Messages: {message_count}")
    print("-" * 40)
    
    # Sequence of operations
    operations = [
        ('P1', 'read', None),
        ('P2', 'read', None),
        ('P1', 'write', 1),
        ('P3', 'read', None),
        ('P2', 'write', 2),
        ('P4', 'read', None),
        ('P1', 'read', None),
    ]

    for i, (proc, op, value) in enumerate(operations):
        print(f"Step {i+1}: {proc} {'writes X = ' + str(value) if op == 'write' else 'reads X'}")
        
        messages_this_step = 0
        message_type = ""
        
        if op == 'read':
            if caches[proc]['state'] in ['I']: # Read Miss
                messages_this_step += 1
                
                # Check other caches
                other_caches_states = [c['state'] for p, c in caches.items() if p != proc]
                
                if 'M' in other_caches_states:
                    message_type = "BusRd (serviced by cache in Modified state)"
                    # Find owner in Modified state
                    owner = next(p for p, c in caches.items() if c['state'] == 'M')
                    read_value = caches[owner]['value']
                    # Owner transitions to Shared and writes to memory
                    caches[owner]['state'] = 'S'
                    memory_value = read_value
                elif 'E' in other_caches_states:
                    message_type = "BusRd (serviced by cache in Exclusive state)"
                    # Find owner in Exclusive state
                    owner = next(p for p, c in caches.items() if c['state'] == 'E')
                    read_value = caches[owner]['value']
                    # Owner transitions to Shared
                    caches[owner]['state'] = 'S'
                elif 'S' in other_caches_states:
                    message_type = "BusRd (serviced by cache in Shared state)"
                    # Find a sharer
                    sharer = next(p for p, c in caches.items() if c['state'] == 'S')
                    read_value = caches[sharer]['value']
                else: # Serviced by memory
                    message_type = "BusRd (serviced by memory)"
                    read_value = memory_value
                    
                # Update current processor's cache
                if 'S' in other_caches_states or 'M' in other_caches_states or 'E' in other_caches_states:
                    caches[proc]['state'] = 'S'
                else:
                    caches[proc]['state'] = 'E'
                caches[proc]['value'] = read_value
            # else Read Hit (no message)
        
        elif op == 'write':
            if caches[proc]['state'] in ['S', 'I']: # Write Miss or needs invalidation
                messages_this_step += 1
                if caches[proc]['state'] == 'S':
                    message_type = "BusUpgr (Invalidate)"
                else: # State is I
                    message_type = "BusRdX (Read Exclusive)"
                
                # Invalidate all other copies
                for p, c in caches.items():
                    if p != proc:
                        if c['state'] in ['S', 'E', 'M']:
                            # If another cache was Modified, it would write back to memory first
                            if c['state'] == 'M':
                                memory_value = c['value']
                            c['state'] = 'I'
                            c['value'] = None
                            
            # else Write Hit on M or E (no bus message)

            # Update current processor's cache to Modified
            caches[proc]['state'] = 'M'
            caches[proc]['value'] = value

        message_count += messages_this_step
        if messages_this_step > 0:
            print(f"Message: {message_type}")
        else:
            print("Message: None (Local cache hit)")

        print(f"Cache States: { {p: c['state'] for p, c in caches.items()} }")
        print(f"Running Total Messages: {message_count}")
        print("-" * 40)

    print(f"The final number of cache coherence messages is {message_count}.")
    
    # Reconstructing the equation for the final output
    equation = "1 + 1 + 1 + 1 + 1 + 1 + 1"
    print(f"Calculation: {equation} = {message_count}")

simulate_mesi_protocol()
<<<7>>>