import math

def calculate_ur_values(start_n, end_n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n)
    for a given range of n.

    The problem specifies the Hamiltonian V(q) = 1/2 * (q^2 - q^n). The period
    integral corresponds to the hyperelliptic curve y^2 = q^n - q^2 + 2E.

    According to a theorem by E. Horozov, for a curve of the form
    y^2 = x^n + a_{n-2}x^{n-2} + ... (where the x^{n-1} term is missing),
    the minimal order of the Picard-Fuchs equation for the period is floor(n/2).
    This applies here.
    """
    
    # We want to find the set {u_r(3), u_r(4), ..., u_r(12)}
    n_values = range(start_n, end_n + 1)
    
    results = []
    for n in n_values:
        # In Python, integer division // is equivalent to the floor function.
        order = n // 2
        results.append(order)
        
    return results

# Set the range for n as specified in the problem
start_n = 3
end_n = 12

# Calculate the list of values
final_values = calculate_ur_values(start_n, end_n)

# Print the final result
print(f"The set {{u_r(3), u_r(4), ..., u_r(12)}} is calculated as:")
print(final_values)