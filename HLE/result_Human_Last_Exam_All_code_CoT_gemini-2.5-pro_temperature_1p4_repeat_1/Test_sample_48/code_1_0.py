try:
    # We use the automata-lib library to perform the necessary conversions.
    # If not installed, you can install it using: pip install automata-lib
    from automata.fa.nfa import NFA
except ImportError:
    print("The 'automata-lib' library is required to run this script.")
    print("Please install it using: 'pip install automata-lib'")
    exit()

# The language L is defined by the regular expression:
# ( (b | c)^* a ( a | b c | c b | b b b )^* (b a | b b a | c a) )^* (b | c)^* a ( a | b c | c b | b b b )^*
# We define this regular expression as a string, removing spaces for the parser.
# The syntax used by the library is standard for regular expressions.
regex_string = "((b|c)*a(a|bc|cb|bbb)*(ba|bba|ca))*(b|c)*a(a|bc|cb|bbb)*"

# Step 1: Create an NFA from the regular expression.
# The NFA.from_regex() class method parses the string and builds the corresponding NFA.
nfa = NFA.from_regex(regex_string)

# Step 2: Convert the NFA to a DFA.
# The .to_dfa() method applies the subset construction algorithm.
dfa = nfa.to_dfa()

# Step 3: Minimize the DFA.
# The .minify() method computes the minimal DFA using Hopcroft's algorithm.
minimal_dfa = dfa.minify()

# Step 4: Count the states of the minimal DFA.
# The number of states is the size of the set of states in the automaton.
num_states = len(minimal_dfa.states)

# Step 5: Print the final result.
# The final result is the number of states in the minimal DFA.
print(f"The number of states in the minimal DFA is: {num_states}")
