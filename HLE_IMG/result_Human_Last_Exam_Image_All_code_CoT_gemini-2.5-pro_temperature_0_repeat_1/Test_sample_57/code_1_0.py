import sys

def fibonacci(k):
    """Calculates the k-th Fibonacci number (F_1=1, F_2=1)."""
    if k <= 0:
        return 0
    if k <= 2:
        return 1
    a, b = 1, 1
    for _ in range(3, k + 1):
        a, b = b, a + b
    return b

def solve_tiling(n):
    """
    Calculates the number of ways to tile the L-shape for a given n
    using the formula 2 * F_{2n-2}.
    """
    if n < 2:
        print(f"The shape is not well-defined for n={n} (must be n >= 2).")
        return

    subscript = 2 * n - 2
    fib_val = fibonacci(subscript)
    result = 2 * fib_val
    
    # Output the equation with the calculated values
    print(f"For n = {n}, the number of ways is given by the equation:")
    print(f"2 * F_({2}*n - {2}) = 2 * F_{subscript} = 2 * {fib_val} = {result}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python your_script_name.py <n>")
        print("Please provide an integer value for n (n >= 2).")
        print("\nThe general formula is: 2 * F_{2n-2}")
    else:
        try:
            n_val = int(sys.argv[1])
            solve_tiling(n_val)
        except ValueError:
            print("Invalid input. Please provide an integer value for n.")
