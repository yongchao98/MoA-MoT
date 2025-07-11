try:
    from automata.fa.nfa import NFA
    from automata.fa.dfa import DFA
except ImportError:
    print("This script requires the 'automata-lib' library.")
    print("Please install it using: pip install automata-lib")
    exit()

# The regular expression provided in the problem.
# Spaces are removed as they are not supported by the library's parser.
regex = "((b|c)*a(a|bc|cb|bbb)*(ba|bba|ca))*(b|c)*a(a|bc|cb|bbb)*"

# Step 1: Construct an NFA from the regular expression.
nfa = NFA.from_regex(regex)

# Step 2: Convert the NFA to a DFA.
dfa = DFA.from_nfa(nfa)

# Step 3: Minimize the DFA.
minimized_dfa = dfa.minify()

# Step 4: Count the number of states in the minimized DFA.
num_states = len(minimized_dfa.states)

# As per the instruction to "output each number in the final equation",
# we print the resulting number in a full sentence.
print(f"The number of states in the minimal DFA is {num_states}.")
