import itertools

def interpreter(program_str):
    """
    A simple interpreter for a primitive recursive RPN-style language.
    This language uses single digits '0'-'9' and operators '+' and '*'.
    Every program is guaranteed to halt. It may fail on malformed programs.
    """
    stack = []
    for char in program_str:
        if char.isdigit():
            stack.append(int(char))
        elif char in ('+', '*'):
            if len(stack) < 2:
                # Not enough operands for an operation, invalid program
                raise ValueError("Invalid program: Not enough operands.")
            op2 = stack.pop()
            op1 = stack.pop()
            if char == '+':
                stack.append(op1 + op2)
            elif char == '*':
                stack.append(op1 * op2)
    
    if len(stack) == 1:
        return stack[0]
    else:
        # Final state should have one number on the stack
        raise ValueError("Invalid program: Stack has more than one value at the end.")

def find_shortest_program(target_n):
    """
    Finds the length of the shortest program that outputs target_n.
    This function implements the brute-force search to compute K(n).
    """
    print(f"Searching for the shortest program to compute n = {target_n}...\n")
    
    # The alphabet of our simple programming language 'P'
    alphabet = '0123456789+*'
    
    # We search programs in increasing order of their length
    for length in itertools.count(1):
        print(f"Checking all {len(alphabet)**length} programs of length {length}...")
        
        # Generate all possible programs of the current length
        possible_programs = itertools.product(alphabet, repeat=length)
        
        for program_tuple in possible_programs:
            program_str = "".join(program_tuple)
            try:
                # Run the program through our interpreter
                result = interpreter(program_str)
                
                # Check if the output matches our target number
                if result == target_n:
                    print("\n--- Found a solution! ---")
                    print(f"The shortest program that outputs {target_n} is: '{program_str}'")
                    print(f"The length of this program is: {length}")
                    # The problem asks to output each number in the final equation.
                    # We will format this as: K(n) = length
                    print(f"K({target_n}) = {length}")
                    return length
            except ValueError:
                # This program was malformed (e.g., '*34'), so we ignore it.
                continue

if __name__ == '__main__':
    # Let's find K(n) for a sample number.
    # A larger number would take significantly longer due to the exponential search space.
    n_to_find = 42
    find_shortest_program(n_to_find)
