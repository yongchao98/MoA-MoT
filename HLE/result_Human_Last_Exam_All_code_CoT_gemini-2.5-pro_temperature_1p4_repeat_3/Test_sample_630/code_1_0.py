import math

def analyze_rates(T):
    """
    Calculates the values of different convergence rates for a given T.
    
    Args:
        T (int): The number of samples.
    """
    
    # Rates from the answer choices
    rate_A = 1 / T
    rate_B = 1 / (T**(2/3))
    rate_C = 1 / (T**(1/2))
    
    # Rate derived from standard convex optimization analysis
    # This corresponds to the minimax optimal rate Theta(log(T) / sqrt(T))
    rate_D_minimax = math.log(T) / math.sqrt(T)
    
    print(f"Analysis for T = {T}:")
    print("-" * 30)
    print(f"Rate A (1/T): \t\t{rate_A:.10f}")
    print(f"Rate B (1/T^(2/3)): \t{rate_B:.10f}")
    print(f"Rate C (1/T^(1/2)): \t{rate_C:.10f}")
    print(f"Derived Rate (log(T)/T^(1/2)): {rate_D_minimax:.10f}")
    print("-" * 30)

    # Illustrating the d-dependent rate from exp-concavity
    # Rate = d * log(T) / T
    # We test for different relationships between d and T
    d_small = 10
    d_medium = int(math.sqrt(T))
    d_large = T
    
    rate_exp_concave_small_d = d_small * math.log(T) / T
    rate_exp_concave_medium_d = d_medium * math.log(T) / T
    rate_exp_concave_large_d = d_large * math.log(T) / T
    
    print("Illustration of the dimension-dependent rate O(d*log(T)/T):")
    # Optimal rate is min(rate_D_minimax, rate_exp_concave)
    
    print(f"For d = {d_small} (small, d < sqrt(T)), optimal rate is ~d*log(T)/T = {min(rate_D_minimax, rate_exp_concave_small_d):.10f}")
    print(f"For d = {d_medium} (medium, d ~ sqrt(T)), optimal rate is ~log(T)/sqrt(T) = {min(rate_D_minimax, rate_exp_concave_medium_d):.10f}")
    print(f"For d = {d_large} (large, d > sqrt(T)), optimal rate is ~log(T)/sqrt(T) = {min(rate_D_minimax, rate_exp_concave_large_d):.10f}")


# Use a large value for T to see the asymptotic behavior
T_samples = 1_000_000
analyze_rates(T_samples)