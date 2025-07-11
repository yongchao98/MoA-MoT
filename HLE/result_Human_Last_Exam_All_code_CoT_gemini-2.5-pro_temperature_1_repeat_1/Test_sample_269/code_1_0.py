def generate_cosmo_loco_output():
    """
    This function simulates the output of the given Scheme code.
    The Scheme code prints 1, then restarts and prints 1 again,
    and finally enters an infinite loop printing 2.
    This script will print the first few terms of this sequence.
    """
    
    # First, 1 is printed.
    print(1, end=' ')
    
    # Then the computation restarts, and 1 is printed again.
    print(1, end=' ')
    
    # Then the code enters an infinite loop, printing 2 repeatedly.
    # We will print it 15 times to show the pattern.
    for _ in range(15):
        print(2, end=' ')
    
    # We add an ellipsis to indicate the sequence continues infinitely.
    print("...")

generate_cosmo_loco_output()