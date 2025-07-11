def solve():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    caches = {
        'P1': {'state': 'I', 'value': None},
        'P2': {'state': 'I', 'value': None},
        'P3': {'state': 'I', 'value': None},
        'P4': {'state': 'I', 'value': None},
    }
    memory_value = 0
    message_count = 0
    processors = ['P1', 'P2', 'P3', 'P4']

    def print_state(op_num, description):
        states = [f"{p}:{caches[p]['state']}" for p in processors]
        print(f"{op_num}. {description}")
        print(f"   States: [{', '.join(states)}], Messages: {message_count}\n")

    # --- Initial State ---
    print("Initial State:")
    print(f"   Memory Value: {memory_value}")
    print_state(0, "All caches are Invalid.")

    # 1. P1 reads X
    # P1 has a read miss, issues a Read Miss message.
    # No other cache has a copy, memory supplies data. P1 -> E
    message_count += 1
    caches['P1']['state'] = 'E'
    caches['P1']['value'] = memory_value
    print_state(1, "P1 reads X (miss). Issues Read Miss. Becomes Exclusive.")

    # 2. P2 reads X
    # P2 has a read miss, issues a Read Miss message.
    # P1 snoops, has X in state E. P1 supplies data. P1: E -> S, P2 -> S.
    message_count += 1
    caches['P1']['state'] = 'S'
    caches['P2']['state'] = 'S'
    caches['P2']['value'] = caches['P1']['value']
    print_state(2, "P2 reads X (miss). Issues Read Miss. P1 provides data. P1->S, P2->S.")

    # 3. P1 writes X = 1
    # P1 has X in state S. To write, it must upgrade. Issues an Invalidate message.
    # P2 invalidates its copy (S -> I). P1 updates and becomes M.
    message_count += 1
    new_value = 1
    caches['P1']['state'] = 'M'
    caches['P1']['value'] = new_value
    caches['P2']['state'] = 'I'
    caches['P2']['value'] = None
    print_state(3, f"P1 writes X={new_value} (hit). Issues Invalidate. P1->M, P2->I.")

    # 4. P3 reads X
    # P3 has a read miss, issues a Read Miss message.
    # P1 snoops, has X in state M. P1 writes back to memory and supplies data.
    # P1: M -> S, P3 -> S.
    message_count += 1
    memory_value = caches['P1']['value']
    caches['P3']['state'] = 'S'
    caches['P3']['value'] = memory_value
    caches['P1']['state'] = 'S'
    print_state(4, "P3 reads X (miss). Issues Read Miss. P1 flushes & provides data. P1->S, P3->S.")

    # 5. P2 writes X = 2
    # P2 has a write miss (state I). Issues RWITM (Read With Intent to Modify).
    # P1 and P3 snoop and invalidate (S -> I). P2 gets data and becomes M.
    message_count += 1
    new_value = 2
    for p in ['P1', 'P3']:
        if caches[p]['state'] != 'I':
            caches[p]['state'] = 'I'
            caches[p]['value'] = None
    caches['P2']['state'] = 'M'
    caches['P2']['value'] = new_value
    print_state(5, f"P2 writes X={new_value} (miss). Issues RWITM. P1 & P3 Invalidate. P2->M.")

    # 6. P4 reads X
    # P4 has a read miss, issues a Read Miss message.
    # P2 snoops (state M). P2 writes back to memory and supplies data.
    # P2: M -> S, P4 -> S.
    message_count += 1
    memory_value = caches['P2']['value']
    caches['P4']['state'] = 'S'
    caches['P4']['value'] = memory_value
    caches['P2']['state'] = 'S'
    print_state(6, "P4 reads X (miss). Issues Read Miss. P2 flushes & provides data. P2->S, P4->S.")

    # 7. P1 reads X
    # P1 has a read miss (state I). Issues a Read Miss message.
    # P2 and P4 are sharing. One of them (or memory) provides the data. P1 becomes S.
    message_count += 1
    caches['P1']['state'] = 'S'
    caches['P1']['value'] = memory_value
    print_state(7, "P1 reads X (miss). Issues Read Miss. Becomes Shared.")

    print(f"Final calculation: 1 + 1 + 1 + 1 + 1 + 1 + 1 = {message_count}")
    print(f"\nTotal cache coherence messages: {message_count}")
    return message_count

final_answer = solve()
print(f"<<<{final_answer}>>>")