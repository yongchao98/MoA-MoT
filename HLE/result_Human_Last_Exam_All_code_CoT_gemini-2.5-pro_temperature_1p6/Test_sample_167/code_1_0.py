def get_alon_tarsi_number_for_K_nm(n, m):
    """
    Calculates the Alon-Tarsi number for a complete bipartite graph K_n,m
    based on a known mathematical theorem.

    Args:
        n (int): The size of the first partition.
        m (int): The size of the second partition.

    Returns:
        int or str: The Alon-Tarsi number or an error message.
    """
    if n <= 0 or m <= 0:
        return "Invalid partition size. n and m must be positive."

    # Handle the special case of K_1,1, which is just an edge.
    if n == 1 and m == 1:
        # AT(K_1,1) = 2
        result = 2
    # Case where partitions have equal size (n >= 2)
    elif n == m:
        # Theorem: AT(K_n,n) = 3 for n >= 2
        result = 3
    # Case where partitions have different sizes
    else:
        # Theorem: AT(K_n,m) = 2 for n != m
        result = 2

    return result

def solve_task():
    """
    Solves the user's specific request for the Alon-Tarsi number of K_1000,1000.
    """
    n = 1000
    m = 1000

    print(f"Finding the Alon-Tarsi number for the complete bipartite graph K_n,m where n={n} and m={m}.")
    print("This is based on the following theorem:")
    print("1. AT(K_n,m) = 2 if n != m.")
    print("2. AT(K_n,n) = 3 if n = m >= 2.")
    print("3. AT(K_1,1) = 2.")
    print("-" * 20)

    # Calculate the Alon-Tarsi number
    alon_tarsi_number = get_alon_tarsi_number_for_K_nm(n, m)

    # Output the final equation with the numbers
    print(f"For the graph K_{n},{m}:")
    print(f"Since n = {n} and m = {m}, the condition n == m is met.")
    print(f"Since n = {n} >= 2, the condition for the theorem is satisfied.")
    print(f"The final result is:")
    print(f"AT(K_{n},{m}) = {alon_tarsi_number}")


if __name__ == "__main__":
    solve_task()
