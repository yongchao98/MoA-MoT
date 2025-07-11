def fibonacci_sequence(n):
    """
    Generates and prints the first n Fibonacci numbers and the final equation.
    """
    if n <= 0:
        print("Please enter a positive integer.")
        return
    
    # First 100 numbers means from F(0) to F(99)
    # The list will hold the sequence.
    fib_numbers = []
    
    a, b = 0, 1
    for _ in range(n):
        fib_numbers.append(a)
        a, b = b, a + b
        
    # Print all the numbers in the sequence
    print("The first 100 Fibonacci numbers are:")
    for num in fib_numbers:
        print(num)
        
    print("\n" + "="*30)
    print("The final equation in the sequence is:")
    
    # For n=100, we have F(0) to F(99).
    # The last number is F(99), which is F(98) + F(97).
    if n >= 3:
        # fib_numbers[-1] is F(99)
        # fib_numbers[-2] is F(98)
        # fib_numbers[-3] is F(97)
        num1 = fib_numbers[-3]
        num2 = fib_numbers[-2]
        result = fib_numbers[-1]
        
        # Output each number in the final equation
        print(f"{num1} + {num2} = {result}")

# Call the function to generate the first 100 Fibonacci numbers.
fibonacci_sequence(100)