# First, ensure the 'automata-lib' library is installed.
# If not, you can install it from your terminal using: pip install automata-lib

try:
    from automata.fa.nfa import NFA
except ImportError:
    print("The required 'automata-lib' library is not found.")
    print("Please install it by running the following command in your terminal:")
    print("pip install automata-lib")
    exit()

# The regular expression for the language L.
# Concatenation is represented by adjacency (no spaces), and '|' for union.
# Parentheses are used for grouping.
regex_str = "(((b|c)*a(a|bc|cb|bbb)*(ba|bba|ca)))*(b|c)*a(a|bc|cb|bbb)*"

try:
    # --- Step 1 & 2: Create an NFA from the regular expression ---
    # The `from_regex` method parses the string and builds an NFA.
    nfa = NFA.from_regex(regex_str)

    # --- Step 3: Convert the NFA to a DFA ---
    # The `to_dfa` method performs the subset construction.
    dfa = nfa.to_dfa()

    # --- Step 4: Minimize the DFA ---
    # The `minify` method returns a new minimal DFA by merging equivalent states.
    minimal_dfa = dfa.minify()

    # --- Step 5: Count the states ---
    # The number of states is the size of the set of states in the minimal DFA.
    num_states = len(minimal_dfa.states)

    # Output the explanation and the result.
    print("Regular Expression:")
    print(f"L = {regex_str}\n")
    print("The number of states in the minimal DFA that recognizes this language can be found by:")
    print("1. Converting the regex to an NFA.")
    print("2. Converting the NFA to a DFA.")
    print("3. Minimizing the DFA.")
    print("\nAfter performing these steps computationally, the result is:\n")
    
    # The user request mentioned a "final equation", which is ambiguous.
    # The most important part is the final number. I will print it clearly.
    print("Final Result:")
    print(f"The number of states in the minimal DFA is {num_states}.")

except Exception as e:
    print(f"An error occurred during the process: {e}")
