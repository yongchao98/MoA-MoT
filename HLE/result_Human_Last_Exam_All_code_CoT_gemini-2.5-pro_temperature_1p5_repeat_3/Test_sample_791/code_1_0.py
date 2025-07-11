import sys

# In Python, the maximum recursion depth is limited.
# We can increase it to show the principle more clearly, but the error is the same.
# Note: A true infinite loop wouldn't error out, it would just hang.
# Python's recursion limit saves us from a true hang.
sys.setrecursionlimit(2000)

# We will use a global counter to track the depth of recursion calls.
recursion_depth = 0

def parse_S():
    """
    A placeholder for parsing the S rule.
    A real parser would try all alternatives: N, ADJ N, N ADJ.
    Let's simulate it trying the 'ADJ N' path first.
    """
    print("Trying to parse S -> 'ADJ N'...")
    # This call will trigger the infinite recursion.
    parse_ADJ()
    # The parser would never reach the part to parse N.

def parse_ADJ():
    """
    A placeholder for parsing the ADJ rule.
    A real parser would try its alternatives: 'red', 'or alike', REC.
    Let's assume it tries the REC alternative, which reveals the problem.
    """
    global recursion_depth
    recursion_depth += 1
    
    print(f"  (Depth {recursion_depth}): Trying to parse ADJ -> REC...")
    if recursion_depth > 10:
      print("  (Depth > 10): Likely infinite recursion, stopping demo.")
      # In a real scenario, this would continue until a stack overflow.
      raise RecursionError("Maximum recursion depth exceeded (simulated)")
      
    parse_REC()

def parse_REC():
    """
    This function demonstrates the left-recursion problem.
    The rule is REC = REC ADJ.
    The function for REC must immediately call the function for REC again.
    """
    print(f"  (Depth {recursion_depth}): Trying to parse REC -> REC ADJ...")
    print(f"  (Depth {recursion_depth}): This requires calling parse_REC() again immediately.")
    
    # This is the left-recursive call that causes the infinite loop.
    parse_REC() # It calls itself without consuming input.
    
    # This code is unreachable
    print("This line will never be printed.")
    parse_ADJ()


def main():
    """
    Main function to run the demonstration and explain the result.
    """
    print("--- Demonstrating the Left-Recursion Problem for RD/PR Parsers ---")
    try:
        # We simulate the start of the parsing process for a sentence.
        # The parser will try to recognize a Subject (S).
        parse_S()
    except RecursionError as e:
        print("\n------------------------------------------------------------")
        print(f"SUCCESS: Caught a RecursionError as expected.")
        print("This error occurs because a function (`parse_REC`) calls itself endlessly")
        print("without making progress through the input string.")
        print("This is the classic failure mode for RD/PR parsers on left-recursive grammars.")
        print("------------------------------------------------------------\n")

    print("--- Evaluating the Answer Choices ---")
    print("\nA. The language contains 'red frogs swim swiftly', and it is not the longest sentence in the language.")
    print("   ANALYSIS: While true for the language itself, it's FALSE in the context of the parsers. The parsers will get stuck in an infinite loop and cannot accept this or any other sentence.")
    
    print("\nB. RD accepts 'snake big jump well.' but PR does not.")
    print("   ANALYSIS: FALSE. The word 'big' is not in the grammar.")
    
    print("\nC. Both RD and PR accept 'jump or alike frog swiftly.' and do not run inifinitely.")
    print("   ANALYSIS: FALSE. The sentence does not match the grammar's S V C structure.")

    print("\nD. Neither RD nor PR can be implemented because the grammar contains left recursion.")
    print("   ANALYSIS: TRUE. 'Implemented' here means creating a working, terminating parser. Due to the infinite loop demonstrated above, this is not possible for these parser types without transforming the grammar, which is disallowed by the problem statement.")
    
    print("\nE. None of the above are correct.")
    print("   ANALYSIS: FALSE, because D is the correct statement.")


if __name__ == "__main__":
    main()