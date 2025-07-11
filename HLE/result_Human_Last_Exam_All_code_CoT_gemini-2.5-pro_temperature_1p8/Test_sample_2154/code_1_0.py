import math

def calculate_ur_sequence(start_n, end_n):
    """
    Calculates the sequence u_r(n) for a given range of n.

    The minimal order u_r(n) of the Picard-Fuchs differential equation
    for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) is given by the formula
    u_r(n) = floor((n-1)/2), according to the work of L. Gavrilov.

    This function computes and prints the sequence for n from start_n to end_n.
    """
    print(f"Calculating the sequence u_r(n) for n = {start_n}, ..., {end_n} using the formula u_r(n) = floor((n-1)/2).\n")
    
    results = []
    for n in range(start_n, end_n + 1):
        # Calculate the intermediate value n-1
        numerator = n - 1
        # Calculate the division
        division_result = numerator / 2
        # Apply the floor function to get the final result
        result = math.floor(division_result)
        
        # Store the result
        results.append(result)
        
        # Print the calculation steps for clarity
        print(f"For n = {n}:")
        print(f"u_r({n}) = floor(({n} - 1) / 2) = floor({numerator} / 2) = floor({division_result}) = {result}\n")
        
    # Print the final list of results
    print("The complete sequence {u_r(3), u_r(4), ..., u_r(12)} is:")
    print(results)


# Define the range for n
start_n = 3
end_n = 12

# Run the calculation and print the results
calculate_ur_sequence(start_n, end_n)
