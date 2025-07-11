# Initial state
a, b, c, d, e, f = 0, 1, 1, 1, 2, 1

# Infinite loop to generate the sequence
while True:
    print a,  # In Python 2, this prints the value of 'a' without a newline

    # Update the state variables
    # This is an integer-based Linear Feedback Shift Register (LFSR)
    # The recurrence relation for 'c' is c[n+4] = c[n] + c[n+1]
    c, d, e, f = d, e, f, c + d
    
    b -= 1
    
    # If the counter 'b' becomes negative, flip the output bit 'a'
    if b < 0:
        a ^= 1
        # and reset 'b' based on the new value of 'c'
        b += 99 * c