import collections

def solve():
    """
    This function defines a sample DFA and a word, and verifies if the word
    synchronizes all reachable states to the unique accept state.

    This code acts as a verifier for the condition for membership in G,
    given a candidate synchronizing word 'w'. The problem of finding 'w'
    is the hard part. If 'w' were guaranteed to have polynomial length,
    a non-deterministic machine could guess 'w' and use this verifier,
    placing the problem in NP.
    """
    
    # Define a Deterministic Finite Automaton (DFA) M
    # M = (Q, Sigma, delta, q0, F)
    # states are represented by integers.
    dfa = {
        'states': {0, 1, 2},
        'alphabet': {'a', 'b'},
        'transitions': {
            (0, 'a'): 1, (0, 'b'): 0,
            (1, 'a'): 1, (1, 'b'): 2,
            (2, 'a'): 2, (2, 'b'): 2,
        },
        'start_state': 0,
        'accept_states': {2}, # Problem specifies one accept state
    }
    
    # Candidate synchronizing word
    word_w = "ab"

    # --- Verification Logic ---
    
    if len(dfa['accept_states']) != 1:
        print("Error: The DFA must have exactly one accept state.")
        return
        
    accept_state = list(dfa['accept_states'])[0]
    
    # 1. Find the set of reachable states R(M) using Breadth-First Search (BFS)
    reachable_states = set()
    q = collections.deque([dfa['start_state']])
    visited = {dfa['start_state']}
    
    while q:
        current_state = q.popleft()
        reachable_states.add(current_state)
        for symbol in dfa['alphabet']:
            next_state = dfa['transitions'].get((current_state, symbol))
            if next_state is not None and next_state not in visited:
                visited.add(next_state)
                q.append(next_state)

    print(f"DFA defined with {len(dfa['states'])} states.")
    print(f"Reachable states found: {sorted(list(reachable_states))}")
    print(f"Unique accept state: {accept_state}")
    print(f"Checking word w = '{word_w}'")

    # 2. For each reachable state, simulate the DFA on word 'w'
    all_sync = True
    for start_node in reachable_states:
        current_node = start_node
        for symbol in word_w:
            current_node = dfa['transitions'].get((current_node, symbol), -1)
            if current_node == -1:
                print(f"Error: Invalid transition from state {current_node} with symbol {symbol}")
                return

        print(f"  - From state {start_node}, word '{word_w}' leads to state {current_node}")
        if current_node != accept_state:
            all_sync = False
            
    # 3. Output the result
    if all_sync:
        print(f"\nResult: True. The word '{word_w}' synchronizes all reachable states to the accept state {accept_state}.")
    else:
        print(f"\nResult: False. The word '{word_w}' does not synchronize all reachable states to the accept state {accept_state}.")

solve()