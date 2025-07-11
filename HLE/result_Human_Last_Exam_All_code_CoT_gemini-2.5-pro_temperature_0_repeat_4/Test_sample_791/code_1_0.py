import sys

# We use a manual depth counter to prevent an actual stack overflow crash.
recursion_depth = 0
max_depth = 15 # A limit to detect the infinite loop.

def parse_rec():
    """
    Simulates the parsing function for the left-recursive REC non-terminal.
    REC = REC ADJ
    """
    global recursion_depth
    recursion_depth += 1
    print(f"Entering parse_rec() at depth {recursion_depth}...")

    if recursion_depth > max_depth:
        print("\n!!! Maximum recursion depth reached. Infinite loop detected! !!!")
        # In a real program, this would likely cause a RecursionError (stack overflow).
        # We exit to stop the simulation.
        sys.exit()

    # This is the direct left-recursive call. The function calls itself
    # without consuming any input, leading to an infinite loop.
    print("   -> From rule 'REC = REC ADJ', calling parse_rec() again.")
    parse_rec()
    # The parser would never reach this point.
    # parse_adj()

def parse_adj():
    """
    Simulates the parsing function for the ADJ non-terminal.
    ADJ = 'red' | 'or alike' | REC
    """
    print("Attempting to parse a non-terminal ADJ...")
    # A real parser would check the input stream here.
    # We simulate the case where the input is neither 'red' nor 'or alike'.
    print("   - Input does not match 'red'.")
    print("   - Input does not match 'or alike'.")

    # The parser now tries the third alternative for ADJ, which is REC.
    print("   - Trying the 'REC' production rule.")
    parse_rec()

# --- Main execution ---
print("Demonstrating the effect of left recursion on a Recursive Descent parser.")
print("The grammar contains the rule 'REC = REC ADJ'.")
print("Let's simulate what happens when the parser tries to recognize an ADJ")
print("that isn't 'red' or 'or alike', forcing it to try the 'REC' rule.")
print("-" * 60)

# This call will trigger the infinite recursion.
parse_adj()
