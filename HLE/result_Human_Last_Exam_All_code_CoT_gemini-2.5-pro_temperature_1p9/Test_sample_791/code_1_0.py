def main():
    """
    This program demonstrates why the left-recursive grammar cannot be handled
    by a standard Recursive Descent or Packrat parser.
    """

    print("Analyzing the grammar for language L:")
    print("L = S V C '.' EOF")
    print("S = N | ADJ N | N ADJ")
    print("N = 'frogs' | 'snakes'")
    print("V = 'jump' | 'swim'")
    print("ADJ = 'red' | 'or alike' | REC")
    print("REC = REC ADJ  <-- This rule is the problem")
    print("C = 'well' | 'swiftly'\n")

    print("The rule 'REC = REC ADJ' is directly left-recursive because the non-terminal")
    print("'REC' appears as the very first symbol on the right side of its own definition.")
    print("\nRecursive Descent (RD) and Packrat (PR) parsers work by creating functions")
    print("for each non-terminal. A function for REC, let's call it parse_REC(),")
    print("would have to call itself as its first action to satisfy the rule.")
    print("This results in infinite recursion, causing the program to crash (stack overflow).")

    print("\nLet's simulate this with a simplified trace:")
    print("---")
    # In a real scenario, this call would be triggered when trying to parse
    # a part of the input as a REC.
    simulate_parse_rec(0)
    print("---\n")

    print("As the simulation shows, any attempt to execute the logic for the REC rule")
    print("leads to a non-terminating loop. Because a parser for L must handle all of its rules,")
    print("and the rule for REC cannot be implemented in a terminating way, a complete and")
    print("correct RD or PR parser for this grammar cannot be implemented.")

    print("\nEvaluating the choices:")
    print("A: While the sentence is in the language, this choice doesn't address the parser implementation failure.")
    print("B & C: The sentences described are not valid according to the grammar.")
    print("D: This correctly identifies the core issue. The left recursion in the grammar makes it impossible to build a working RD or PR parser without modification.")
    print("\nThe correct statement is D.")

def simulate_parse_rec(depth):
    """A function to simulate the infinite recursion of parsing a left-recursive rule."""
    max_depth = 4
    indent = "  " * depth
    print(f"{indent}Calling parse_REC()...")

    if depth >= max_depth:
        print(f"{indent}  -> Reached simulation limit! In a real parser, this would continue until a stack overflow error.")
        return

    # According to the rule 'REC = REC ADJ', the first step is to call parse_REC() again.
    print(f"{indent}  -> First step of rule 'REC = REC ADJ' is to parse REC.")
    simulate_parse_rec(depth + 1)


if __name__ == '__main__':
    main()
<<<D>>>