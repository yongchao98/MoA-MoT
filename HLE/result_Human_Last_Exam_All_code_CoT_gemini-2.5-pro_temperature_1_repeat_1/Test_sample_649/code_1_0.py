import math

def calculate_minimax_lower_bound(n, delta, N, kl_divergences, phi_val_at_delta_half):
    """
    Calculates the minimax lower bound based on Tsybakov's method.

    The formula is: R*_n >= (Phi(delta/2)/2) * (1 - sqrt((n * K) / 2))
    where K is the average KL divergence: K = (1/N) * sum(D_KL(P_j || P_0)).

    Args:
        n (int): The number of i.i.d. datapoints.
        delta (float): The minimum separation between the null parameter and alternative parameters.
        N (int): The number of alternative hypotheses.
        kl_divergences (list[float]): A list of KL divergences [D_KL(P_1||P_0), ..., D_KL(P_N||P_0)].
        phi_val_at_delta_half (float): The value of the loss function Phi at delta/2.
                                       For example, if Phi(x) = x^2, this would be (delta/2)**2.
    """
    if len(kl_divergences) != N:
        raise ValueError("The length of kl_divergences list must be equal to N.")
    
    # Calculate K, the average KL divergence
    avg_kl = sum(kl_divergences) / N
    
    # Calculate the term inside the square root
    sqrt_inner_term = (n * avg_kl) / 2
    
    # The bound is non-trivial only if this term is less than 1.
    if sqrt_inner_term >= 1:
        print("The condition for a non-trivial bound (n * K / 2 < 1) is not met.")
        print(f"Calculated value for n * K / 2: {sqrt_inner_term:.4f}")
        lower_bound = 0.0
    else:
        # Calculate the main parenthesis term
        parenthesis_term = 1 - math.sqrt(sqrt_inner_term)
        
        # Calculate the final lower bound
        lower_bound = (phi_val_at_delta_half / 2) * parenthesis_term

    # Print the equation with the final numbers
    print("Derived Lower Bound Formula:")
    print("R*_n >= (Phi(delta/2) / 2) * (1 - sqrt(n * K / 2))")
    print("\nSubstituting the given values:")
    print(f"R*_n >= ({phi_val_at_delta_half:.4f} / 2) * (1 - sqrt({n} * {avg_kl:.4f} / 2))")
    print(f"R*_n >= {phi_val_at_delta_half/2:.4f} * (1 - sqrt({sqrt_inner_term:.4f}))")
    if sqrt_inner_term < 1:
        print(f"R*_n >= {phi_val_at_delta_half/2:.4f} * (1 - {math.sqrt(sqrt_inner_term):.4f})")
        print(f"R*_n >= {phi_val_at_delta_half/2:.4f} * {1 - math.sqrt(sqrt_inner_term):.4f}")
    
    final_result = max(0, lower_bound)
    print(f"\nFinal tightest lower bound proved: {final_result:.4f}")
    
    # The final answer in the required format
    print(f"\n<<<{final_result:.4f}>>>")


if __name__ == '__main__':
    # --- User-configurable parameters ---
    
    # Number of i.i.d. datapoints
    n_samples = 100
    
    # Number of alternative hypotheses
    N_hypotheses = 10
    
    # Minimum separation of parameters
    delta_separation = 0.5
    
    # The value of Phi(delta/2). For example, if Phi(x) = x^2, this would be (0.5/2)^2 = 0.0625.
    # If Phi(x) = x, this would be 0.25. Let's use Phi(x)=x for this example.
    phi_value = delta_separation / 2
    
    # A list of KL-divergences D_KL(P_j || P_0) for j = 1 to N.
    # For this example, let's assume they are all the same small value.
    # The value must be small enough so that n*K/2 < 1.
    # Let's choose a value like 0.0003 for each.
    kl_divergences_list = [0.0003] * N_hypotheses
    
    # --- Calculation ---
    calculate_minimax_lower_bound(
        n=n_samples,
        delta=delta_separation,
        N=N_hypotheses,
        kl_divergences=kl_divergences_list,
        phi_val_at_delta_half=phi_value
    )