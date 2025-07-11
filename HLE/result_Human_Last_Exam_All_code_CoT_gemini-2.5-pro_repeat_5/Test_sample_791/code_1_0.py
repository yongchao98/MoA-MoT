import sys

# Set a higher recursion limit to show more steps, but not run forever.
# A real infinite recursion would crash without this.
sys.setrecursionlimit(20)

def demonstrate_left_recursion():
    """
    This function simulates the behavior of a Recursive Descent parser
    on the given grammar to demonstrate the problem with left recursion.
    """
    print("--- Parser Simulation Start ---")
    print("Simulating a parser trying to parse a rule involving 'REC'.")
    print("The problematic rules are 'ADJ = ... | REC' and 'REC = REC ADJ'.\n")

    # We use a call_depth tracker to visualize the recursion.
    call_depth = 0

    # Define the parsing functions based on the grammar.

    def parse_rec(depth):
        """Simulates parsing the REC non-terminal."""
        indent = "  " * depth
        print(f"{indent}Enter parse_rec() at depth {depth}.")

        # The grammar rule is REC = REC ADJ.
        # The parser immediately tries to parse REC again.
        print(f"{indent}Rule 'REC = REC ADJ' makes function call itself immediately.")
        
        # This is the infinitely recursive call.
        # It happens before any input is consumed or any other rule is checked.
        parse_rec(depth + 1)

    def parse_adj(depth):
        """Simulates parsing the ADJ non-terminal."""
        indent = "  " * depth
        print(f"{indent}Enter parse_adj() at depth {depth}.")
        
        # A real parser would try to match terminals first.
        # Let's assume input doesn't match 'red' or 'or alike'.
        print(f"{indent}Assume input is not 'red' or 'or alike'.")
        print(f"{indent}Parser will now try the 'REC' alternative in 'ADJ = ... | REC'.")
        
        # Now, call the function for REC.
        parse_rec(depth + 1)

    try:
        # Start the simulation by trying to parse an ADJ.
        parse_adj(call_depth)
    except RecursionError:
        print("\n--- Simulation Halted ---")
        print("A 'RecursionError' was caught.")
        print("This demonstrates that the left-recursive rule 'REC = REC ADJ' causes the parser to enter an infinite loop of function calls.")
        print("Standard Recursive Descent and Packrat parsers cannot handle this.")
        print("\nConclusion: Neither RD nor PR can be correctly implemented for this grammar due to left recursion.")

demonstrate_left_recursion()