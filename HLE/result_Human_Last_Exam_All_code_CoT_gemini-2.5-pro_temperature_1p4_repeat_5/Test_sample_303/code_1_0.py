def fibonacci(n):
    """
    Calculates the nth Fibonacci number using recursion.
    This function demonstrates the correct structure with base cases.
    """
    # Base cases: The recursion stops here.
    # This is what the "Caramel Sauce" recipe is missing.
    if n <= 0:
        return 0
    elif n == 1:
        return 1
    # Recursive step: The function calls itself with smaller inputs.
    else:
        return fibonacci(n-1) + fibonacci(n-2)

def main():
    """
    Main function to print the first 20 Fibonacci numbers.
    The "equation" is the sequence of numbers itself.
    """
    # The original "recipe" aimed for 100 numbers, but we will show the first 20.
    num_terms = 20
    print(f"The first {num_terms} Fibonacci numbers are:")
    
    # We will print each number of the sequence.
    # This corresponds to "output each number in the final equation!"
    fib_sequence = []
    for i in range(num_terms):
        fib_sequence.append(str(fibonacci(i)))
    
    print(", ".join(fib_sequence))

if __name__ == "__main__":
    main()
