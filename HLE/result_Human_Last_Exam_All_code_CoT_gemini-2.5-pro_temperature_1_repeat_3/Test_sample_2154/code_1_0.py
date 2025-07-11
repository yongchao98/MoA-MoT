import math

def calculate_u_r_sequence():
    """
    Calculates the sequence u_r(n) for n from 3 to 12.
    
    u_r(n) is the minimal order of the Picard-Fuchs differential equation for the
    period of the Hamiltonian V(q) = 1/2 * (q^2 - q^n).
    
    The formula depends on the parity of n:
    - If n is odd, u_r(n) = n - 1.
    - If n is even, u_r(n) = floor(n/4 - 1) + 1.
    """
    
    u_r_values = []
    for n in range(3, 13):
        if n % 2 != 0:
            # Odd n
            order = n - 1
        else:
            # Even n
            order = math.floor(n / 4 - 1) + 1
        u_r_values.append(order)
    
    # The final question asks for the set of values.
    # We will print the list representing this set.
    print(u_r_values)

calculate_u_r_sequence()