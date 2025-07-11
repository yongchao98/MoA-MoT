# The 'greenery' library is required for this script.
# You can install it via pip: pip install greenery
try:
    from greenery import lego
except ImportError:
    print("Please install the 'greenery' library to run this code:")
    print("pip install greenery")
    exit()

# The language L is defined by the regular expression:
# ( (b | c)^* a ( a | b c | c b | b b b )^* (b a | b b a | c a) )^* (b | c)^* a ( a | b c | c b | b b b )^*
regex = r"(((b|c)*a(a|bc|cb|bbb)*(ba|bba|ca))*(b|c)*a(a|bc|cb|bbb)*)"

# Parse the regular expression into a pattern object.
pattern = lego.parse(regex)

# Convert the pattern to a finite-state machine (FSM).
# This results in a non-deterministic automaton (NFA).
fsm = pattern.to_fsm()

# The .reduce() method converts the NFA to a DFA and minimizes it.
minimal_dfa = fsm.reduce()

# The number of states in the minimal DFA.
num_states = len(minimal_dfa.states)

# The question asks to output the number in the final equation.
# The number of states is the final answer.
print(num_states)