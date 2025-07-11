def simulate_mesi_coherence():
    """
    Simulates the MESI protocol for a sequence of operations and counts coherence messages.
    """
    # States for variable X in caches of P1, P2, P3, P4
    states = ['I', 'I', 'I', 'I']
    # Value of X in memory and caches
    memory_value = 0
    cache_values = [0] * 4
    
    total_messages = 0
    message_counts_per_step = []

    def print_state(step_description, messages_this_step, explanation):
        nonlocal total_messages
        total_messages += messages_this_step
        message_counts_per_step.append(messages_this_step)
        
        print(f"--- {step_description} ---")
        print("\n".join(explanation))
        print(f"\nMessages for this step: {messages_this_step}")
        print(f"New States: P1: {states[0]}, P2: {states[1]}, P3: {states[2]}, P4: {states[3]}")
        
        equation_str = " + ".join(map(str, message_counts_per_step))
        print(f"Total Messages so far: {equation_str} = {total_messages}\n")

    print("Initial State:")
    print(f"P1: {states[0]}, P2: {states[1]}, P3: {states[2]}, P4: {states[3]}")
    print(f"Value of X in Memory: {memory_value}")
    print(f"Total Messages: {total_messages}\n")

    # Step 1: P1 reads X
    states[0] = 'E'
    cache_values[0] = memory_value
    print_state("Step 1: P1 reads X", 1, [
        "- P1 has a Read Miss, as its state is Invalid.",
        "- P1 issues a 'Read Miss' message on the bus.",
        "- Memory provides the data (0). Since P1 is the only cache holding the data, its state becomes Exclusive."
    ])
    
    # Step 2: P2 reads X
    states[0] = 'S'
    states[1] = 'S'
    cache_values[1] = cache_values[0]
    print_state("Step 2: P2 reads X", 1, [
        "- P2 has a Read Miss.",
        "- P2 issues a 'Read Miss' message.",
        "- P1 snoops the bus, sees the request for X, and provides the data via a cache-to-cache transfer.",
        "- Both P1 and P2 now share the data, so their states become Shared."
    ])
    
    # Step 3: P1 writes X = 1
    states[0] = 'M'
    states[1] = 'I'
    cache_values[0] = 1
    print_state("Step 3: P1 writes X = 1", 1, [
        "- P1 must write, but its state is Shared. It cannot write until other copies are invalidated.",
        "- P1 broadcasts an 'Invalidate' message on the bus.",
        "- P2 snoops and invalidates its copy (S -> I).",
        "- P1 updates its value and its state becomes Modified."
    ])

    # Step 4: P3 reads X
    states[0] = 'S'
    states[2] = 'S'
    cache_values[2] = cache_values[0]
    memory_value = cache_values[0] # Write-back to memory on M -> S transition
    print_state("Step 4: P3 reads X", 1, [
        "- P3 has a Read Miss and issues a 'Read Miss' message.",
        "- P1 snoops, sees it has a Modified copy, and supplies the data (1) to P3.",
        "- P1 also writes the data back to main memory.",
        "- P1 and P3 now share the data, so their states become Shared."
    ])

    # Step 5: P2 writes X = 2
    states[1] = 'M'
    states[0] = 'I'
    states[2] = 'I'
    cache_values[1] = 2
    print_state("Step 5: P2 writes X = 2", 1, [
        "- P2 has a Write Miss (state is Invalid).",
        "- P2 broadcasts a 'Read-for-Ownership' (RFO) message to get the data and invalidate others.",
        "- P1 and P3 snoop the RFO, invalidate their copies (S -> I), and one supplies the data.",
        "- P2 receives the data, updates its value to 2, and its state becomes Modified."
    ])

    # Step 6: P4 reads X
    states[1] = 'S'
    states[3] = 'S'
    cache_values[3] = cache_values[1]
    memory_value = cache_values[1] # Write-back to memory
    print_state("Step 6: P4 reads X", 1, [
        "- P4 has a Read Miss and issues a 'Read Miss' message.",
        "- P2 snoops, provides its Modified data (2), and writes back to memory.",
        "- P2 and P4 now share the data, so their states become Shared."
    ])
    
    # Step 7: P1 reads X
    states[0] = 'S'
    cache_values[0] = cache_values[1]
    print_state("Step 7: P1 reads X", 1, [
        "- P1 has a Read Miss and issues a 'Read Miss' message.",
        "- P2 or P4, which have Shared copies, will provide the data via a cache-to-cache transfer.",
        "- P1's state becomes Shared."
    ])

    print("--- Final Calculation ---")
    equation_str = " (P1 Read) + ".join(map(str, message_counts_per_step[0:2]))
    equation_str += " (P1 Write) + " + str(message_counts_per_step[2])
    equation_str += " (P3 Read) + " + str(message_counts_per_step[3])
    equation_str += " (P2 Write) + " + str(message_counts_per_step[4])
    equation_str += " (P4 Read) + " + str(message_counts_per_step[5])
    equation_str += " (P1 Read)"

    final_equation = f"1 (P1 Read) + 1 (P2 Read) + 1 (P1 Write) + 1 (P3 Read) + 1 (P2 Write) + 1 (P4 Read) + 1 (P1 Read) = {total_messages}"
    print(final_equation)
    print(f"\nTotal cache coherence messages exchanged: {total_messages}")


if __name__ == "__main__":
    simulate_mesi_coherence()