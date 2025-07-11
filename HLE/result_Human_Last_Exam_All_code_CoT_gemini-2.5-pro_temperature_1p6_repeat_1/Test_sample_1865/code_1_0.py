def solve_mesi_coherence():
    """
    Calculates the number of cache coherence messages for a given sequence of operations
    in a 4-processor system using the MESI protocol.
    """
    
    # Each element in the list contains an explanation for a step and the number of messages.
    steps = [
        ("P1 reads X", "P1 has a read miss and issues a BusRd message. Since no other cache has X, memory supplies the data and P1's cache state becomes Exclusive (E).", 1),
        ("P2 reads X", "P2 has a read miss and issues a BusRd message. P1 snoops the bus, supplies the data, and its state changes from Exclusive to Shared (S). P2's state becomes S.", 1),
        ("P1 writes X = 1", "P1 holds X in state S and needs to write. It sends an Invalidate message to gain exclusive ownership. P2 invalidates its copy. P1's state becomes Modified (M).", 1),
        ("P3 reads X", "P3 has a read miss and issues a BusRd message. P1 snoops and finds its copy is Modified. P1 must perform a Write-Back to memory and supply the data to P3. This requires two bus transactions (BusRd and Write-Back). P1's state becomes S.", 2),
        ("P2 writes X = 2", "P2 has a write miss (its copy is Invalid) and issues a BusRdX message. This single message fetches the data and invalidates the Shared copies in P1 and P3. P2's state becomes Modified (M).", 1),
        ("P4 reads X", "P4 has a read miss and issues a BusRd message. P2 snoops and finds its copy is Modified. Similar to step 4, P2 performs a Write-Back and supplies the data. This requires two bus transactions. P2's state becomes S.", 2),
        ("P1 reads X", "P1 has a read miss and issues a BusRd message. P2 and P4 hold clean copies in the Shared state, so one of them supplies the data. No write-back is needed.", 1)
    ]
    
    total_messages = 0
    message_counts = []
    
    print("Analyzing the sequence of operations:\n")
    
    for i, (operation, explanation, msg_count) in enumerate(steps):
        print(f"Step {i+1}: {operation}")
        print(f"  - Action: {explanation}")
        print(f"  - Coherence Messages: {msg_count}\n")
        total_messages += msg_count
        message_counts.append(str(msg_count))

    print("-" * 30)
    print("Final Calculation")
    print("-" * 30)
    equation = " + ".join(message_counts)
    print(f"The total number of messages is the sum of messages from each step:")
    print(f"{equation} = {total_messages}")

solve_mesi_coherence()
<<<9>>>