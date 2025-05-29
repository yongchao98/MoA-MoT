def complete_parentheses(s):
    # Count opening and closing brackets
    open_count = s.count('(')
    close_count = s.count(')')
    
    # Calculate how many closing brackets we need
    needed_close = open_count - close_count
    
    # Add the needed closing brackets
    result = s + ')' * needed_close
    
    # Verify if it's valid
    def is_valid(seq):
        count = 0
        for char in seq:
            if char == '(':
                count += 1
            else:
                count -= 1
            if count < 0:
                return False
        return count == 0
    
    if is_valid(result):
        print(f"Complete sequence: {result}")
        return True
    else:
        print("Cannot create valid sequence")
        return False

# Test the input
input_str = "( ( ( ) )"
complete_parentheses(input_str)