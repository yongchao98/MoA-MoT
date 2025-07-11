def solve_mesi_simulation():
    """
    Simulates the MESI protocol for a given sequence of operations and
    counts the number of cache coherence messages.
    """
    # caches[processor] = [state, value]
    caches = {
        'P1': ['I', None], 'P2': ['I', None],
        'P3': ['I', None], 'P4': ['I', None]
    }
    memory_value = 0
    message_count = 0

    # A helper to print the current state
    def print_state():
        state_str = ", ".join([f"P{i+1}: {caches[f'P{i+1}'][0]}" for i in range(4)])
        print(f"  Current States: {state_str}")
        print("-" * 40)

    print("Initial State:")
    print(f"  Memory X = {memory_value}")
    print_state()

    # --- Operation 1: P1 reads X ---
    message_count += 1
    print(f"1. P1 reads X (Value={memory_value})")
    print(f"   - P1 issues a 'Read Miss'. (Messages so far: {message_count})")
    caches['P1'] = ['E', memory_value]
    print(f"   - P1 becomes Exclusive as it's the only cache with X.")
    print_state()

    # --- Operation 2: P2 reads X ---
    message_count += 1
    print(f"2. P2 reads X (Value={caches['P1'][1]})")
    print(f"   - P2 issues a 'Read Miss'. (Messages so far: {message_count})")
    caches['P1'][0] = 'S' # P1 state changes from E to S
    caches['P2'] = ['S', caches['P1'][1]] # P2 receives data
    print(f"   - P1 provides data. P1->S, P2->S.")
    print_state()

    # --- Operation 3: P1 writes X = 1 ---
    message_count += 1
    print(f"3. P1 writes X = 1")
    print(f"   - P1 sends an 'Invalidate' to upgrade S->M. (Messages so far: {message_count})")
    caches['P1'] = ['M', 1]
    caches['P2'][0] = 'I' # P2 is invalidated
    print(f"   - P2 invalidates its copy. P1 becomes Modified.")
    print_state()

    # --- Operation 4: P3 reads X ---
    message_count += 1
    print(f"4. P3 reads X (Value={caches['P1'][1]})")
    print(f"   - P3 issues a 'Read Miss'. (Messages so far: {message_count})")
    memory_value = caches['P1'][1] # P1 writes back to memory
    caches['P1'][0] = 'S' # P1 state changes from M to S
    caches['P3'] = ['S', memory_value] # P3 receives data
    print(f"   - P1 provides data (and writes to memory). P1->S, P3->S.")
    print_state()

    # --- Operation 5: P2 writes X = 2 ---
    message_count += 1
    print(f"5. P2 writes X = 2")
    print(f"   - P2 issues a 'Read For Ownership' (RFO). (Messages so far: {message_count})")
    caches['P1'][0] = 'I' # P1 is invalidated
    caches['P3'][0] = 'I' # P3 is invalidated
    caches['P2'] = ['M', 2]
    print(f"   - P1 and P3 are invalidated. P2 becomes Modified.")
    print_state()

    # --- Operation 6: P4 reads X ---
    message_count += 1
    print(f"6. P4 reads X (Value={caches['P2'][1]})")
    print(f"   - P4 issues a 'Read Miss'. (Messages so far: {message_count})")
    memory_value = caches['P2'][1] # P2 writes back to memory
    caches['P2'][0] = 'S' # P2 state changes from M to S
    caches['P4'] = ['S', memory_value] # P4 receives data
    print(f"   - P2 provides data (and writes to memory). P2->S, P4->S.")
    print_state()

    # --- Operation 7: P1 reads X ---
    message_count += 1
    print(f"7. P1 reads X (Value={caches['P2'][1]})")
    print(f"   - P1 issues a 'Read Miss'. (Messages so far: {message_count})")
    caches['P1'] = ['S', memory_value] # P1 receives data from another cache
    print(f"   - P2 or P4 provides data. P1 becomes Shared.")
    print_state()

    print("\nFinal Calculation:")
    final_equation = "1 + 1 + 1 + 1 + 1 + 1 + 1"
    print(f"Total cache coherence messages = {final_equation} = {message_count}")

solve_mesi_simulation()