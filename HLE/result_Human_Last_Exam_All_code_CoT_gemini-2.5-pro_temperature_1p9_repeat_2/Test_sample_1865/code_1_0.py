def calculate_mesi_messages():
    """
    Simulates a sequence of memory operations in a 4-processor system
    using the MESI protocol and counts the cache coherence messages.
    """
    # caches[processor] = state of the cache line for variable X
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    total_messages = 0
    message_counts_per_step = []
    
    # Helper function to print the state of the simulation
    def print_step(step, processor, operation, description, messages_this_step, new_states):
        nonlocal total_messages
        print(f"--- Step {step}: {processor} {operation} ---")
        for line in description:
            print(f"  - {line}")
        
        # Update states
        for p, s in new_states.items():
            caches[p] = s
        
        print(f"  - New States: P1:{caches['P1']}, P2:{caches['P2']}, P3:{caches['P3']}, P4:{caches['P4']}")
        print(f"  - Messages This Step: {messages_this_step}")
        total_messages += messages_this_step
        message_counts_per_step.append(messages_this_step)
        print("-" * 35 + "\n")

    print("Initial State: All caches are Invalid (I).\n")

    # 1. P1 reads X
    print_step(
        step=1, processor="P1", operation="reads X",
        description=[
            "P1 has a read miss and sends a 'Read Miss' message.",
            "Memory supplies the data since no other cache has it.",
            "P1 becomes the exclusive owner of the data.",
        ],
        messages_this_step=1,
        new_states={'P1': 'E'}
    )
    
    # 2. P2 reads X
    print_step(
        step=2, processor="P2", operation="reads X",
        description=[
            "P2 has a read miss and sends a 'Read Miss' message.",
            "P1 snoops the bus, supplies the data (cache-to-cache transfer).",
            "Both P1 and P2 now share the data.",
        ],
        messages_this_step=1,
        new_states={'P1': 'S', 'P2': 'S'}
    )

    # 3. P1 writes X = 1
    print_step(
        step=3, processor="P1", operation="writes X = 1",
        description=[
            "P1 needs to write but only has a Shared copy.",
            "P1 sends an 'Invalidate' message to gain exclusive ownership.",
            "P2 invalidates its copy. P1's copy becomes Modified.",
        ],
        messages_this_step=1,
        new_states={'P1': 'M', 'P2': 'I'}
    )

    # 4. P3 reads X
    print_step(
        step=4, processor="P3", operation="reads X",
        description=[
            "P3 has a read miss and sends a 'Read Miss' message.",
            "P1 snoops, supplies the modified data, and writes it back to memory.",
            "P1 and P3 now share the updated data.",
        ],
        messages_this_step=1,
        new_states={'P1': 'S', 'P3': 'S'}
    )

    # 5. P2 writes X = 2
    print_step(
        step=5, processor="P2", operation="writes X = 2",
        description=[
            "P2 has a write miss (its copy is Invalid).",
            "P2 sends a 'Read-For-Ownership' (RFO) message to get the data and invalidate others.",
            "P1 and P3 invalidate their copies. P2's copy becomes Modified.",
        ],
        messages_this_step=1,
        new_states={'P1': 'I', 'P3': 'I', 'P2': 'M'}
    )

    # 6. P4 reads X
    print_step(
        step=6, processor="P4", operation="reads X",
        description=[
            "P4 has a read miss and sends a 'Read Miss' message.",
            "P2 snoops, supplies the modified data, and writes it back to memory.",
            "P2 and P4 now share the updated data.",
        ],
        messages_this_step=1,
        new_states={'P2': 'S', 'P4': 'S'}
    )
    
    # 7. P1 reads X
    print_step(
        step=7, processor="P1", operation="reads X",
        description=[
            "P1 has a read miss and sends a 'Read Miss' message.",
            "A sharer (P2 or P4) supplies the data via cache-to-cache transfer.",
            "P1's state becomes Shared.",
        ],
        messages_this_step=1,
        new_states={'P1': 'S'}
    )

    # Final Calculation
    print("--- Total Messages Calculation ---")
    equation_str = " + ".join(map(str, message_counts_per_step))
    print(f"The total number of messages is the sum of messages from each step:")
    print(f"{equation_str} = {total_messages}")


# Run the simulation
if __name__ == "__main__":
    calculate_mesi_messages()
    print("\n<<<7>>>")