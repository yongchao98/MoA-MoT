import math

def calculate_picard_fuchs_orders():
    """
    This function calculates the minimal order of the Picard-Fuchs differential
    equation, u_r(n), for the given Hamiltonian for n from 3 to 12.

    The formula for the order is u_r(n) = 2 * g, where g is the genus of the
    associated hyperelliptic curve. The genus g is given by floor((n-1)/2).
    """
    
    print("Calculating the sequence {u_r(3), u_r(4), ..., u_r(12)}:")
    
    results = []
    # Iterate through the given range of n
    for n in range(3, 13):
        # Calculate u_r(n) using the formula 2 * floor((n-1)/2).
        # In Python, integer division `//` achieves the floor operation.
        order = 2 * ((n - 1) // 2)
        results.append(order)
        # The prompt requires printing each number of the "final equation".
        # We interpret this as printing each term of the sequence.
        print(f"u_r({n}) = {order}")

if __name__ == "__main__":
    calculate_picard_fuchs_orders()
