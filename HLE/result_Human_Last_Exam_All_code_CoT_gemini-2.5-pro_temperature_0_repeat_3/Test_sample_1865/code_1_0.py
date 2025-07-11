def solve_mesi_messages():
    """
    Simulates the MESI protocol for a sequence of operations and counts the messages.
    """
    # Initial states
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    processors = list(caches.keys())
    total_messages = 0
    message_log = []

    # Helper function to print state
    def print_state(step_description, step_messages):
        nonlocal total_messages
        total_messages += step_messages
        message_log.append(step_messages)
        print(f"{step_description}")
        print(f"  - States: {', '.join([f'{p}({s})' for p, s in caches.items()])}")
        print(f"  - Messages in this step: {step_messages}")
        print("-" * 20)

    # --- Simulation Start ---
    print("Initial State: All caches are Invalid (I).\n" + "-" * 20)

    # 1. P1 reads X
    # Read miss, no other cache has it. BusRd + Data from Memory.
    caches['P1'] = 'E'
    print_state("1. P1 reads X", 2)

    # 2. P2 reads X
    # Read miss, P1 has it in E. BusRd + Data from P1. P1 -> S.
    caches['P1'] = 'S'
    caches['P2'] = 'S'
    print_state("2. P2 reads X", 2)

    # 3. P1 writes X = 1
    # Write hit, but state is S. BusUpgr/Invalidate. P2 -> I. P1 -> M.
    caches['P1'] = 'M'
    caches['P2'] = 'I'
    print_state("3. P1 writes X = 1", 1)

    # 4. P3 reads X
    # Read miss, P1 has it in M. BusRd + Data from P1 (and WB to memory). P1 -> S.
    caches['P1'] = 'S'
    caches['P3'] = 'S'
    print_state("4. P3 reads X", 2)

    # 5. P2 writes X = 2
    # Write miss (state I). BusRdX + Data from P1/P3/Memory. P1,P3 -> I. P2 -> M.
    caches['P1'] = 'I'
    caches['P3'] = 'I'
    caches['P2'] = 'M'
    print_state("5. P2 writes X = 2", 2)

    # 6. P4 reads X
    # Read miss, P2 has it in M. BusRd + Data from P2 (and WB to memory). P2 -> S.
    caches['P2'] = 'S'
    caches['P4'] = 'S'
    print_state("6. P4 reads X", 2)

    # 7. P1 reads X
    # Read miss, P2/P4 have it in S. BusRd + Data from another cache.
    caches['P1'] = 'S'
    print_state("7. P1 reads X", 2)

    # Final Calculation
    equation = " + ".join(map(str, message_log))
    print(f"Total messages exchanged: {equation} = {total_messages}")

solve_mesi_messages()
<<<13>>>