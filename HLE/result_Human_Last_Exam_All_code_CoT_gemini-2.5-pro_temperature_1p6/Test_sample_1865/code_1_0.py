def simulate_mesi_protocol():
    """
    Simulates the MESI protocol for a given sequence of operations and counts coherence messages.
    Each bus transaction (BusRd, BusUpgr, BusRdX) is counted as one message.
    """
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    cache_values = {'P1': 0, 'P2': 0, 'P3': 0, 'P4': 0}
    memory_value = 0
    message_count = 0
    message_breakdown = []

    operations = [
        ('P1', 'read', None),
        ('P2', 'read', None),
        ('P1', 'write', 1),
        ('P3', 'read', None),
        ('P2', 'write', 2),
        ('P4', 'read', None),
        ('P1', 'read', None),
    ]

    print("--- MESI Protocol Simulation ---")
    print(f"Initial State: Caches={caches}, Memory={memory_value}\n")

    for i, (proc, op, value) in enumerate(operations):
        print(f"Step {i+1}: {proc} {'writes' if op == 'write' else 'reads'} X" + (f" with value {value}" if value is not None else ""))
        
        current_state = caches[proc]
        message_sent = False
        message_type = ""

        if op == 'read':
            if current_state == 'I': # Read Miss
                message_sent = True
                message_type = "Read Miss (BusRd)"
                
                # Check other caches
                source = 'memory'
                data_value = memory_value
                for p, state in caches.items():
                    if state in ['E', 'M', 'S']:
                        source = p
                
                # If another cache has it in M state
                if caches.get(source) == 'M':
                    print(f"  - {proc} has a read miss. Issues {message_type}.")
                    print(f"  - {source} (in M) snoops, writes back to memory, and supplies data.")
                    memory_value = cache_values[source]
                    data_value = cache_values[source]
                    caches[source] = 'S'
                # If another cache has it in E or S state
                elif caches.get(source) in ['E', 'S']:
                    print(f"  - {proc} has a read miss. Issues {message_type}.")
                    print(f"  - {source} (in {caches[source]}) snoops and supplies data.")
                    data_value = cache_values[source]
                    if caches[source] == 'E':
                        caches[source] = 'S'
                # If only memory has it
                else:
                    print(f"  - {proc} has a read miss. Issues {message_type}.")
                    print("  - Memory supplies the data.")

                caches[proc] = 'S'
                # Special case: if no other cache had the block, state becomes Exclusive
                is_shared = any(p != proc and s != 'I' for p, s in caches.items())
                if not is_shared:
                     caches[proc] = 'E'

                cache_values[proc] = data_value

            else: # Read Hit
                print(f"  - {proc} has a read hit. No message sent.")
                
        elif op == 'write':
            if current_state in ['M', 'E']: # Write Hit (no message)
                print(f"  - {proc} has a write hit in state {current_state}. No message sent.")
                caches[proc] = 'M'
                cache_values[proc] = value

            elif current_state == 'S': # Write Hit (but needs to invalidate)
                message_sent = True
                message_type = "Invalidate (BusUpgr)"
                print(f"  - {proc} writes to a Shared line. Issues {message_type}.")
                for p, state in caches.items():
                    if p != proc and state == 'S':
                        caches[p] = 'I'
                        print(f"  - {p} invalidates its copy.")
                caches[proc] = 'M'
                cache_values[proc] = value
                
            elif current_state == 'I': # Write Miss
                message_sent = True
                message_type = "Read With Intent to Modify (BusRdX)"
                print(f"  - {proc} has a write miss. Issues {message_type}.")
                
                # Invalidate all other copies
                for p, state in caches.items():
                    if p != proc and state != 'I':
                        # If a Modified copy exists, it must write back first
                        if state == 'M':
                           memory_value = cache_values[p]
                           print(f"  - {p} (in M) snoops, writes back to memory, then invalidates.")
                        else:
                           print(f"  - {p} (in {state}) snoops and invalidates its copy.")
                        caches[p] = 'I'

                caches[proc] = 'M'
                cache_values[proc] = value

        if message_sent:
            message_count += 1
            message_breakdown.append(1)

        print(f"  Result -> States: {caches}, Memory Value: {memory_value}")
        print("-" * 25)

    print("\n--- Final Calculation ---")
    calculation_str = " + ".join(map(str, message_breakdown))
    print(f"Total messages = {calculation_str} = {message_count}")


if __name__ == '__main__':
    simulate_mesi_protocol()