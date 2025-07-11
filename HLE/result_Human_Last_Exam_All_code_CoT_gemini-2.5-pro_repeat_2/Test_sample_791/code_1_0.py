import sys

def parse_rec(input_string):
    """
    A mock parsing function for the left-recursive rule: REC = REC ADJ
    """
    # To parse a REC, the first thing it does is try to parse a REC.
    # This function calls itself without consuming any input, leading to an infinite loop.
    print("-> Calling parse_rec...")
    return parse_rec(input_string)

def demonstrate_left_recursion():
    """
    This program demonstrates why a standard Recursive Descent or Packrat parser
    fails on a grammar with left recursion.

    The grammar contains the rule `REC = REC ADJ`. A parser attempting to handle
    this rule would create a function, let's call it `parse_rec`. To satisfy
    the rule, the `parse_rec` function's first action would be to call itself,
    `parse_rec`, without consuming any input. This results in an infinite
    chain of self-calls, which in most programming languages, including Python,
    causes a stack overflow or a recursion depth error.

    The problem explicitly forbids grammar transformations or algorithm modifications
    that would fix this. Therefore, a functional, terminating parser for the
    entire grammar cannot be built using these methods.

    Statement D, "Neither RD nor PR can be implemented because the grammar
    contains left recursion," accurately describes this situation. The following
    code will attempt the infinite recursion and catch the resulting error.
    """
    # Set a low recursion limit to make the demonstration quick and clear.
    # The default is often 1000.
    original_limit = sys.getrecursionlimit()
    low_limit = 50
    sys.setrecursionlimit(low_limit)
    
    print(demonstrate_left_recursion.__doc__)
    print(f"\nSetting a low recursion limit ({low_limit}) to demonstrate the issue quickly.")
    print("Now, simulating the parser's attempt to parse the REC rule:")

    try:
        # This call will trigger the infinite recursion.
        parse_rec("some input string")
    except RecursionError as e:
        print("\n-------------------------------------------------------------")
        print(f"SUCCESS (in demonstrating the failure): Caught an expected RecursionError.")
        print(f"Error message: {e}")
        print("This error proves that the parser entered an infinite loop, as predicted.")
        print("This confirms that statement D is the correct analysis.")
        print("-------------------------------------------------------------")
    finally:
        # Restore the original recursion limit.
        sys.setrecursionlimit(original_limit)

if __name__ == "__main__":
    demonstrate_left_recursion()

<<<D>>>