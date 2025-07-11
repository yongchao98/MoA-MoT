def solve_mesi_coherence():
    """
    Calculates the number of cache coherence messages for a given sequence of operations
    in a 4-processor system using the MESI protocol.
    """
    # Cache states for variable X in each processor P1, P2, P3, P4
    # States: 'M' (Modified), 'E' (Exclusive), 'S' (Shared), 'I' (Invalid)
    cache_states = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    memory_value = 0
    total_messages = 0
    
    # A list to store the breakdown of messages for the final equation
    message_breakdown = []

    print("Initial State:")
    print(f"  Cache States: {cache_states}")
    print(f"  Memory Value of X: {memory_value}")
    print("-" * 40)

    # 1. P1 reads X
    step = 1
    messages_step_1 = 0
    print(f"Step {step}: P1 reads X")
    # P1 has a read miss. It issues a Read Miss (BusRd) to the bus.
    # No other cache has the data, so memory responds.
    # P1's cache state for X becomes Exclusive (E).
    messages_step_1 += 1  # P1: BusRd
    total_messages += messages_step_1
    message_breakdown.append(messages_step_1)
    cache_states['P1'] = 'E'
    print(f"  - P1 issues a Read Miss (BusRd). Memory provides the data.")
    print(f"  - P1's state for X changes from I to E.")
    print(f"  - Messages this step: {messages_step_1} (BusRd)")
    print(f"  - Current States: {cache_states}")
    print("-" * 40)

    # 2. P2 reads X
    step = 2
    messages_step_2 = 0
    print(f"Step {step}: P2 reads X")
    # P2 has a read miss and issues a BusRd.
    # P1 snoops the bus and sees the request. Its state is E.
    # P1 provides the data to P2 via a cache-to-cache transfer (Flush).
    # Both P1 and P2 change their state to Shared (S).
    messages_step_2 += 1  # P2: BusRd
    messages_step_2 += 1  # P1: Flush (data response)
    total_messages += messages_step_2
    message_breakdown.append(messages_step_2)
    cache_states['P1'] = 'S'
    cache_states['P2'] = 'S'
    print(f"  - P2 issues a Read Miss (BusRd).")
    print(f"  - P1 (in state E) provides the data to P2 via a Flush.")
    print(f"  - P1 state changes E -> S. P2 state changes I -> S.")
    print(f"  - Messages this step: {messages_step_2} (BusRd + Flush)")
    print(f"  - Current States: {cache_states}")
    print("-" * 40)

    # 3. P1 writes X = 1
    step = 3
    messages_step_3 = 0
    print(f"Step {step}: P1 writes X = 1")
    # P1's state is S. To write, it must invalidate other shared copies.
    # P1 issues an Invalidate message on the bus.
    # P2 receives the invalidate and changes its state to I.
    # P1's state changes to Modified (M).
    messages_step_3 += 1  # P1: Invalidate
    total_messages += messages_step_3
    message_breakdown.append(messages_step_3)
    cache_states['P1'] = 'M'
    cache_states['P2'] = 'I'
    print(f"  - P1 (in state S) issues an Invalidate message to gain exclusive ownership.")
    print(f"  - P2 invalidates its copy (S -> I). P1's state becomes M (S -> M).")
    print(f"  - Messages this step: {messages_step_3} (Invalidate)")
    print(f"  - Current States: {cache_states}")
    print("-" * 40)

    # 4. P3 reads X
    step = 4
    messages_step_4 = 0
    print(f"Step {step}: P3 reads X")
    # P3 has a read miss and issues a BusRd.
    # P1 snoops and sees its state is M. It must provide the data.
    # P1 flushes the data to P3 and also writes it back to memory.
    # P1 and P3 states become S.
    messages_step_4 += 1  # P3: BusRd
    messages_step_4 += 1  # P1: Flush
    messages_step_4 += 1  # P1: Write-Back to memory
    total_messages += messages_step_4
    message_breakdown.append(messages_step_4)
    cache_states['P1'] = 'S'
    cache_states['P3'] = 'S'
    memory_value = 1
    print(f"  - P3 issues a Read Miss (BusRd).")
    print(f"  - P1 (in state M) flushes data to P3 and writes it back to memory.")
    print(f"  - P1 state changes M -> S. P3 state changes I -> S.")
    print(f"  - Messages this step: {messages_step_4} (BusRd + Flush + Write-Back)")
    print(f"  - Current States: {cache_states}")
    print("-" * 40)

    # 5. P2 writes X = 2
    step = 5
    messages_step_5 = 0
    print(f"Step {step}: P2 writes X = 2")
    # P2's state is I (write miss). It needs data and ownership.
    # P2 issues a Read For Ownership (BusRdX).
    # P1 and P3 snoop. One of them (e.g., P1) provides the data via Flush.
    # Both P1 and P3 invalidate their copies (S -> I).
    # P2's state becomes M.
    messages_step_5 += 1  # P2: BusRdX
    messages_step_5 += 1  # P1 or P3: Flush
    total_messages += messages_step_5
    message_breakdown.append(messages_step_5)
    cache_states['P2'] = 'M'
    cache_states['P1'] = 'I'
    cache_states['P3'] = 'I'
    print(f"  - P2 issues a Read For Ownership (BusRdX).")
    print(f"  - A sharer (P1 or P3) provides the data via Flush.")
    print(f"  - P1 and P3 invalidate their copies (S -> I). P2's state becomes M.")
    print(f"  - Messages this step: {messages_step_5} (BusRdX + Flush)")
    print(f"  - Current States: {cache_states}")
    print("-" * 40)

    # 6. P4 reads X
    step = 6
    messages_step_6 = 0
    print(f"Step {step}: P4 reads X")
    # P4 has a read miss and issues a BusRd.
    # P2 snoops (state M) and provides the data.
    # P2 flushes data to P4 and writes it back to memory.
    # P2 and P4 states become S.
    messages_step_6 += 1  # P4: BusRd
    messages_step_6 += 1  # P2: Flush
    messages_step_6 += 1  # P2: Write-Back
    total_messages += messages_step_6
    message_breakdown.append(messages_step_6)
    cache_states['P2'] = 'S'
    cache_states['P4'] = 'S'
    memory_value = 2
    print(f"  - P4 issues a Read Miss (BusRd).")
    print(f"  - P2 (in state M) flushes data to P4 and writes it back to memory.")
    print(f"  - P2 state changes M -> S. P4 state changes I -> S.")
    print(f"  - Messages this step: {messages_step_6} (BusRd + Flush + Write-Back)")
    print(f"  - Current States: {cache_states}")
    print("-" * 40)

    # 7. P1 reads X
    step = 7
    messages_step_7 = 0
    print(f"Step {step}: P1 reads X")
    # P1 has a read miss and issues a BusRd.
    # P2 and P4 are in state S. One of them (e.g., P2) provides the data via Flush.
    # P1's state becomes S. P2 and P4 remain S.
    messages_step_7 += 1  # P1: BusRd
    messages_step_7 += 1  # P2 or P4: Flush
    total_messages += messages_step_7
    message_breakdown.append(messages_step_7)
    cache_states['P1'] = 'S'
    print(f"  - P1 issues a Read Miss (BusRd).")
    print(f"  - A sharer (P2 or P4) provides the data via Flush.")
    print(f"  - P1's state changes I -> S.")
    print(f"  - Messages this step: {messages_step_7} (BusRd + Flush)")
    print(f"  - Final States: {cache_states}")
    print("-" * 40)

    # Final Calculation
    equation = " + ".join(map(str, message_breakdown))
    print("Total Messages Calculation:")
    print(f"  Sum of messages from each step: {equation} = {total_messages}")
    
    return total_messages

if __name__ == '__main__':
    final_answer = solve_mesi_coherence()
    print(f"\nTotal cache coherence messages exchanged: {final_answer}")
    print(f"<<<{final_answer}>>>")
