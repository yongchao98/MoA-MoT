import sys

# Python's recursion limit is high, but we can lower it for a quick demonstration.
# A standard RD parser loop would be truly infinite.
sys.setrecursionlimit(20)

def demonstrate_parser_failure():
    """
    This function demonstrates the failure of a recursive descent parser
    when it encounters a left-recursive grammar rule.

    The grammar in question contains:
    ADJ = 'red' | 'or alike' | REC
    REC = REC ADJ
    """
    print("--- Parser Demonstration ---")
    print("Analyzing grammar with left recursion: 'REC = REC ADJ'")

    def parse_REC():
        """Simulates the parser function for the left-recursive REC rule."""
        print("-> Calling parse_REC()... This will immediately call itself again.")
        # This recursive call causes an infinite loop.
        parse_REC()

    def parse_ADJ(token):
        """
        Simulates parsing an adjective. If the token is not a known
        terminal, the parser will try the 'REC' alternative.
        """
        print(f"\nAttempting to parse '{token}' as an ADJ.")
        if token == "red":
            print("   Match found: 'red'.")
            return
        
        print("   No terminal match. Trying the 'REC' production rule.")
        parse_REC()

    try:
        # We simulate the parser trying to parse an input like "green frogs..."
        # It would first attempt to parse "green" as an adjective.
        parse_ADJ("green")
    except RecursionError:
        print("\n--- RESULT ---")
        print("PYTHON ERROR: RecursionError: maximum recursion depth exceeded.")
        print("\nThis error proves the concept:")
        print("The function 'parse_REC' called itself without consuming input, creating an infinite loop.")
        print("This is the exact reason why standard Recursive Descent (and Packrat) parsers cannot handle left-recursive grammars.")
        print("\nTherefore, a correct and terminating parser for the given grammar 'L' cannot be implemented using these methods.")

if __name__ == "__main__":
    demonstrate_parser_failure()
