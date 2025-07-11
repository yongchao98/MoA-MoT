try:
    # Ensure the 'automata-lib' library is installed.
    # You can install it using pip: pip install automata-lib
    from automata.fa.nfa import NFA
except ImportError:
    print("The 'automata-lib' library is not installed.")
    print("Please install it by running: pip install automata-lib")
    exit()

def solve_dfa_states():
    """
    This function calculates the number of states in the minimal DFA
    for the given regular expression.
    """
    # The language L is defined by the regular expression:
    # ( (b | c)^* a ( a | b c | c b | b b b )^* (b a | b b a | c a) )^* (b | c)^* a ( a | b c | c b | b b b )^*
    # We represent this as a string that the library can parse.
    # The library uses standard regex syntax: '|' for union, '*' for Kleene star,
    # and juxtaposition for concatenation.
    regex_string = "((b|c)*a(a|bc|cb|bbb)*(ba|bba|ca))*(b|c)*a(a|bc|cb|bbb)*"

    print(f"Regular Expression: {regex_string}")

    # 1. Create an NFA from the regular expression.
    # The library will automatically infer the alphabet {a, b, c}.
    nfa = NFA.from_regex(regex_string)

    # 2. Convert the NFA to a DFA.
    dfa = nfa.to_dfa()

    # 3. Minimize the DFA.
    min_dfa = dfa.minify()

    # 4. Count the number of states in the minimal DFA.
    num_states = len(min_dfa.states)

    print(f"The number of states in the minimal DFA is: {num_states}")
    
    # Although not requested, for completeness, we can show the states and transitions
    # print("\n--- Minimal DFA Details ---")
    # print("States:", min_dfa.states)
    # print("Initial State:", min_dfa.initial_state)
    # print("Final States:", min_dfa.final_states)
    # print("Transitions:")
    # for from_state, transitions in min_dfa.transitions.items():
    #     for symbol, to_state in transitions.items():
    #         print(f"  {from_state} --{symbol}--> {to_state}")
    # print("-------------------------\n")


if __name__ == '__main__':
    solve_dfa_states()
