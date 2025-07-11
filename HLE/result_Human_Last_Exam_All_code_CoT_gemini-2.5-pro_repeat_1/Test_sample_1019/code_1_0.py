import math

def calculate_hamiltonicity_threshold():
    """
    This script calculates the d-threshold for Hamiltonicity for the union of a dense graph and a random graph.

    The problem considers a graph H_n with minimum degree d >= n/2 - eta,
    where 1/2 <= eta <= n/64. The threshold is the probability p such that
    H_n union G(n, p) is Hamiltonian asymptotically almost surely.

    The theoretical threshold is p = Theta((eta/n)^2). We will use the formula p = (eta/n)^2.
    """
    
    # We choose example values for n and eta that fit the problem's constraints.
    # n: number of vertices.
    n = 2000
    
    # eta: parameter for the minimum degree d = n/2 - eta.
    # Constraint: 1/2 <= eta <= n/64. For n=2000, n/64 = 31.25.
    # We choose a value in this range.
    eta = 20
    
    # The minimum degree d for the host graph H_n.
    d = n / 2 - eta
    
    print(f"Calculating the d-threshold for Hamiltonicity for n = {n} and eta = {eta}.")
    print(f"The minimum degree of the host graph is d >= {d}.")
    print("-" * 30)

    # The formula for the threshold probability is p = (eta / n)^2.
    # The prompt requires outputting each number in the final equation.
    
    base_numerator = eta
    base_denominator = n
    power = 2
    
    # Calculate the threshold p
    p = (base_numerator / base_denominator) ** power
    
    print("The threshold probability 'p' is calculated using the formula:")
    print("p = (eta / n) ^ power")
    print("\nThe numbers in this equation for our example are:")
    print(f"eta = {base_numerator}")
    print(f"n = {base_denominator}")
    print(f"power = {power}")
    
    print("\nSubstituting these values:")
    print(f"p = ({base_numerator} / {base_denominator}) ^ {power}")
    
    print(f"\nThe calculated d-threshold is: p = {p}")

# Execute the calculation
calculate_hamiltonicity_threshold()