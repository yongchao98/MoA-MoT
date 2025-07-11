import sys
# It's recommended to install the library 'automata-lib' to run this code.
# You can install it by running: pip install automata-lib
try:
    from automata.fa.nfa import NFA
    from automata.fa.dfa import DFA
except ImportError:
    print("The 'automata-lib' library is not installed.")
    print("Please install it using: pip install automata-lib")
    sys.exit(1)

def solve():
    """
    This function builds the minimal DFA for the given regular expression
    and prints the number of states.
    """
    # Step 1: Define basic NFAs for each symbol in the alphabet {a, b, c}
    # Note: In older versions of the library, NFA.from_symbol was used.
    # In newer versions (>=5.0.0), it's NFA.from_string.
    try:
        nfa_a = NFA.from_string('a')
        nfa_b = NFA.from_string('b')
        nfa_c = NFA.from_string('c')
    except AttributeError:
        # Fallback for older library versions
        nfa_a = NFA.from_symbol('a')
        nfa_b = NFA.from_symbol('b')
        nfa_c = NFA.from_symbol('c')

    # Define a helper for the star operation to handle API changes
    def kleene_star(nfa):
        if hasattr(nfa, 'star'): # For newer versions (>=5.0.0)
            return nfa.star()
        else: # For older versions
            return nfa.kleene_star()

    # Step 2: Build NFAs for the sub-expressions X, Y, Z

    # Build NFA for X = (b | c)^*
    nfa_b_or_c = nfa_b | nfa_c
    nfa_X = kleene_star(nfa_b_or_c)

    # Build NFA for Y = (a | bc | cb | bbb)^*
    nfa_bc = nfa_b.concatenate(nfa_c)
    nfa_cb = nfa_c.concatenate(nfa_b)
    nfa_bbb = nfa_b.concatenate(nfa_b).concatenate(nfa_b)
    nfa_Y_base = nfa_a | nfa_bc | nfa_cb | nfa_bbb
    nfa_Y = kleene_star(nfa_Y_base)

    # Build NFA for Z = (ba | bba | ca)
    nfa_ba = nfa_b.concatenate(nfa_a)
    nfa_bba = nfa_b.concatenate(nfa_b).concatenate(nfa_a)
    nfa_ca = nfa_c.concatenate(nfa_a)
    nfa_Z = nfa_ba | nfa_bba | nfa_ca

    # Step 3: Build the NFA for the full language L = (X a Y Z)^* X a Y

    # Build NFA for the inner block (X a Y Z)
    nfa_Xa = nfa_X.concatenate(nfa_a)
    nfa_XaY = nfa_Xa.concatenate(nfa_Y)
    nfa_XaYZ = nfa_XaY.concatenate(nfa_Z)

    # Apply Kleene star to the inner block
    nfa_XaYZ_star = kleene_star(nfa_XaYZ)

    # Build NFA for the final block (X a Y)
    # We can reuse nfa_XaY from above
    
    # Concatenate the parts to get the final NFA for L
    nfa_L = nfa_XaYZ_star.concatenate(nfa_XaY)

    # Step 4: Convert the NFA for L to a DFA
    dfa_L = DFA.from_nfa(nfa_L)

    # Step 5: Minimize the DFA
    min_dfa_L = dfa_L.minify()

    # Step 6: Count the states and print the result
    num_states = len(min_dfa_L.states)
    print(num_states)

solve()