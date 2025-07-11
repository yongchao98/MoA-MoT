def generate_fibonacci_sequence(n):
    """
    Generates the first n Fibonacci numbers and prints the sequence.
    It also prints the equation for the final number in the sequence.
    """
    if n <= 0:
        print("Please enter a positive integer for the number of terms.")
        return
    
    fib_sequence = []
    a, b = 0, 1
    
    for i in range(n):
        fib_sequence.append(a)
        # Store the previous two numbers for the final equation printout
        if i == n - 2:
            f_n_minus_2 = a
        if i == n - 1:
            f_n_minus_1 = a
            
        a, b = b, a + b
        
    # Per the instruction to "output each number in the final equation",
    # we show the calculation for the final term of the sequence.
    # The 100th term is F(99), which is F(97) + F(98).
    # In our sequence, these are the 98th and 99th elements (index 97 and 98).
    final_num = fib_sequence[n-2] + fib_sequence[n-1]
    
    print("The final equation to calculate the 100th Fibonacci number (F_99) is:")
    print(f"{fib_sequence[n-2]} + {fib_sequence[n-1]} = {final_num}")
    print("-" * 30)

    print(f"The first {n} Fibonacci numbers are:")
    # Print the full list of numbers
    for number in fib_sequence:
        print(number)

# The recipe intended to calculate the first 100 Fibonacci numbers.
generate_fibonacci_sequence(100)