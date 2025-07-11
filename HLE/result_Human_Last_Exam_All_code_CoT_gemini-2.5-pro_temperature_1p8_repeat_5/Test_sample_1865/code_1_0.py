def solve_mesi_messages():
    """
    Calculates the number of MESI coherence messages for a given sequence of operations.
    """
    # Messages for each step in the sequence
    messages_step_1 = 1  # P1 reads X: Read miss, BusRd message
    messages_step_2 = 1  # P2 reads X: Read miss, BusRd message
    messages_step_3 = 1  # P1 writes X: Write hit on Shared state, Invalidate message
    messages_step_4 = 1  # P3 reads X: Read miss, BusRd message
    messages_step_5 = 1  # P2 writes X: Write miss, BusRdX message
    messages_step_6 = 1  # P4 reads X: Read miss, BusRd message
    messages_step_7 = 1  # P1 reads X: Read miss, BusRd message

    # Sum of all messages
    total_messages = (
        messages_step_1 +
        messages_step_2 +
        messages_step_3 +
        messages_step_4 +
        messages_step_5 +
        messages_step_6 +
        messages_step_7
    )

    # Print the equation
    print(f"Total Messages = {messages_step_1} (P1 Rd) + {messages_step_2} (P2 Rd) + {messages_step_3} (P1 Wr) + {messages_step_4} (P3 Rd) + {messages_step_5} (P2 Wr) + {messages_step_6} (P4 Rd) + {messages_step_7} (P1 Rd)")
    print(f"Total Messages = {total_messages}")

solve_mesi_messages()