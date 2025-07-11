import math

def calculate_m_n(n):
    """
    Calculates the minimum number of rewiring operations m(n) needed to
    create a core of high-degree vertices, as described by option J.
    
    Args:
        n (int): The total number of vertices in the network.
        
    Returns:
        tuple: A tuple containing m(n), option_b_val (n), and option_h_val (n/6).
    """
    if n <= 1:
        return 0, 0, 0

    # The maximum degree allowed for any vertex.
    max_degree = math.ceil(math.log2(n))
    
    # The number of vertices in the core, modeled as Θ(log n).
    # For this calculation, we'll use a constant factor of 1.
    num_hubs = math.log2(n)
    
    # The initial degree of each vertex.
    initial_degree = 6
    
    # We only proceed if the target degree is higher than the initial degree.
    if max_degree <= initial_degree:
        return 0, n, n / 6
        
    # Each rewiring operation consists of one edge removal and one edge addition.
    # The total increase in degree for the hubs is the number of hubs multiplied
    # by the degree increase per hub.
    total_degree_increase = num_hubs * (max_degree - initial_degree)
    
    # The number of rewiring operations, m(n), is half of the total degree
    # increase required, as each new edge adds 2 to the sum of degrees.
    m_n = 0.5 * total_degree_increase
    
    return m_n, n, n/6

def main():
    """
    Main function to execute the analysis for a sample network size.
    """
    n = 1_000_000 # A sample large n
    
    max_degree = math.ceil(math.log2(n))
    num_hubs = math.log2(n)
    initial_degree = 6

    m_n, n_val, n_over_6_val = calculate_m_n(n)
    
    print(f"Analyzing the transformation for a network of n = {n:,} vertices.")
    print("-" * 50)
    print("Hypothesis: The transformation requires creating a dense core of high-degree vertices (Option J).")
    print("Let's calculate the minimum cost 'm(n)' for this.")
    print("\n1. Core parameters:")
    print(f"   - Number of hub vertices (k = log_2(n)): {num_hubs:.2f}")
    print(f"   - Target degree for hubs (d = ceil(log_2(n))): {max_degree}")
    print(f"   - Initial degree (k_0): {initial_degree}")

    print("\n2. Cost calculation `m(n) = 0.5 * k * (d - k_0)`:")
    # Here we output each number in the final equation, as requested.
    print(f"   m({n:,}) = 0.5 * {num_hubs:.2f} * ({max_degree} - {initial_degree})")
    print(f"   m({n:,}) = {m_n:.2f}")

    print("\n3. Comparison with other options:")
    print(f"   - m(n) from our calculation = {m_n:,.2f}")
    print(f"   - Value for Θ(n) scaling (Option B) = {n_val:,.2f}")
    print(f"   - Lower bound from Option H (n/6) = {n_over_6_val:,.2f}")
    
    print("\nConclusion:")
    print("The calculated minimum cost m(n) is Θ((log n)^2), which is much smaller than Θ(n).")
    print("This suggests that options B and H are incorrect.")
    print("Therefore, the most accurate statement describes the necessary structure of the final graph.")
    print("Option J provides the best description of this required structure.")


if __name__ == "__main__":
    main()
