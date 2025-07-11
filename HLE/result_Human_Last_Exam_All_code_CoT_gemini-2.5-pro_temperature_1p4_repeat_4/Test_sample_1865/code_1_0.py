def solve_mesi_coherence():
    """
    Simulates the MESI protocol for a given sequence of operations
    and calculates the total number of cache coherence messages.
    """
    # Cache states for variable X in each processor
    # I: Invalid, S: Shared, E: Exclusive, M: Modified
    caches = {'P1': 'I', 'P2': 'I', 'P3': 'I', 'P4': 'I'}
    
    # List to store the number of messages for each step
    message_log = []

    def log_step(operation, new_states, messages_this_step, explanation):
        """Helper function to print the state and update counts."""
        print(f"--- {operation} ---")
        print(f"Explanation: {explanation}")
        
        # Update cache states
        caches.update(new_states)
        
        # Log message count for the final equation
        message_log.append(messages_this_step)
        
        print(f"Messages this step: {messages_this_step}")
        print(f"Cache States: P1={caches['P1']}, P2={caches['P2']}, P3={caches['P3']}, P4={caches['P4']}\n")

    print("Initial State: All caches are Invalid (I).\n")

    # 1. P1 reads X
    # This is a Read Miss. P1 sends a Read request on the bus. Since no other cache has X,
    # memory responds. P1 loads X and its cache line state becomes Exclusive (E).
    log_step("1. P1 reads X", 
             {'P1': 'E'}, 
             1, 
             "P1 issues a Read Miss. No other cache has the data. P1 loads from memory and enters the Exclusive (E) state.")

    # 2. P2 reads X
    # This is a Read Miss for P2. P2 sends a Read request. P1 snoops the bus and sees it holds
    # X in state E. P1 provides the data to P2. Both caches transition to Shared (S).
    log_step("2. P2 reads X", 
             {'P1': 'S', 'P2': 'S'}, 
             1,
             "P2 issues a Read Miss. P1 snoops and provides the data. P1 and P2 caches both enter the Shared (S) state.")

    # 3. P1 writes X = 1
    # This is a Write Hit for P1, but its state is S. P1 must send an Invalidate message on the
    # bus to gain exclusive ownership. P2 receives the invalidation and moves its copy to Invalid (I).
    # P1 writes the new value and its state becomes Modified (M).
    log_step("3. P1 writes X = 1", 
             {'P1': 'M', 'P2': 'I'}, 
             1,
             "P1 has the data in Shared state. It broadcasts an Invalidate message. P2 invalidates its copy. P1's state becomes Modified (M).")

    # 4. P3 reads X
    # This is a Read Miss for P3. P3 sends a Read request. P1 snoops, holding X in state M.
    # P1 flushes the data to memory and also sends it to P3. Both P1 and P3 now have a clean copy,
    # so their states become Shared (S).
    log_step("4. P3 reads X",
             {'P1': 'S', 'P3': 'S'}, 
             1,
             "P3 issues a Read Miss. P1 has the data in Modified state. P1 provides the data and both enter the Shared (S) state.")

    # 5. P2 writes X = 2
    # This is a Write Miss for P2 (state is I). P2 sends a Read-Exclusive (or RWITM) request on the bus.
    # This single message serves to both get the data and invalidate other copies. P1 and P3 receive this
    # and invalidate their copies (S -> I). P2 writes the new value and its state becomes Modified (M).
    log_step("5. P2 writes X = 2", 
             {'P1': 'I', 'P2': 'M', 'P3': 'I'}, 
             1,
             "P2 issues a Read-Exclusive request. P1 and P3 invalidate their copies. P2 writes the data and enters the Modified (M) state.")

    # 6. P4 reads X
    # This is a Read Miss for P4. P4 sends a Read request. P2 snoops, holding X in state M.
    # P2 flushes the data to memory and sends it to P4. Both P2 and P4 now hold a clean copy,
    # so their states become Shared (S).
    log_step("6. P4 reads X", 
             {'P2': 'S', 'P4': 'S'}, 
             1,
             "P4 issues a Read Miss. P2 has the data in Modified state. P2 provides the data and both enter the Shared (S) state.")

    # 7. P1 reads X
    # This is a Read Miss for P1 (state is I). P1 sends a Read request. P2 and P4 both have the data in
    # state S. One of the caches (e.g., P2) responds with the data (cache-to-cache transfer). P1's
    # state becomes Shared (S). P2 and P4 remain in S.
    log_step("7. P1 reads X", 
             {'P1': 'S'}, 
             1,
             "P1 issues a Read Miss. P2 (or P4) has the data in Shared state and provides it. P1 enters the Shared (S) state.")

    # Final calculation
    total_messages = sum(message_log)
    equation = " + ".join(map(str, message_log))
    print("--- Final Calculation ---")
    print(f"Total Messages = {equation} = {total_messages}")

solve_mesi_coherence()
<<<7>>>