import math

def calculate_u_r_n():
    """
    Calculates and prints the minimal order of the Picard-Fuchs differential equation
    u_r(n) for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) for n from 3 to 12.
    """
    
    print("The formula for the minimal order is u_r(n) = 2 * floor((n-1)/2).")
    print("This can be simplified to u_r(n) = n-1 for odd n, and u_r(n) = n-2 for even n.")
    print("-" * 30)
    
    results = []
    
    # Loop for n from 3 to 12
    for n in range(3, 13):
        # In Python, integer division `//` on positive numbers is equivalent to floor.
        order = 2 * ((n - 1) // 2)
        results.append(order)
        
        # Print the detailed calculation for each n, as requested.
        # This illustrates "output each number in the final equation".
        print(f"For n = {n}:")
        print(f"u_r({n}) = 2 * floor(({n}-1)/2) = 2 * floor({(n-1)/2:.1f}) = 2 * {math.floor((n-1)/2)} = {order}")
        print() # Add a newline for better readability

    print("-" * 30)
    print("The complete sequence {u_r(3), u_r(4), ..., u_r(12)} is:")
    print(results)

# Execute the function to display the output
calculate_u_r_n()