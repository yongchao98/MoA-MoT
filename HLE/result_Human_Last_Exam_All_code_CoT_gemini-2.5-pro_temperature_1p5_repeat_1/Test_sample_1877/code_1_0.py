import math

def solve_dynamical_problem():
    """
    This script calculates T(n,p) for the logistic map and derives the final answer.
    """
    
    # --- Introduction and Plan ---
    print("This script solves for the number of iterations T(n,p) to distinguish a periodic from a chaotic orbit in the logistic map.")
    print("The method uses the concept of the Lyapunov exponent to determine the 'predictability horizon' of a chaotic system.")
    print("The formula is T(p) ≈ (p * ln(10)) / λ, where p is the precision and λ is the Lyapunov exponent.")
    print("-" * 50)

    # --- Step 1: Define parameters from the problem ---
    n = 3  # We are analyzing the region near the period-3 window.
    p = 7  # Precision for Era B (7 significant digits).
    
    # We choose a parameter 'r' in the chaotic region near the period-3 window.
    # The period-3 window begins around r ≈ 3.8284. A standard choice for chaos in this vicinity is r=3.83.
    r = 3.83

    print("Step 1: Set up parameters based on the problem description.")
    print(f"  - n = {n} (period to distinguish)")
    print(f"  - p = {p} (precision for Era B)")
    print(f"  - r = {r} (chosen parameter for chaos near the period-3 window)")
    print("-" * 50)
    
    # --- Step 2: Calculate the Lyapunov exponent (λ) ---
    print("Step 2: Numerically calculate the Lyapunov exponent (λ) for r = 3.83.")
    
    x = 0.5  # Standard initial condition
    transient_iterations = 2000  # Iterations to let the system settle onto the attractor
    calc_iterations = 100000     # Iterations to average over for a stable lambda value

    # Discard transient iterations to ensure the orbit is on the attractor
    for _ in range(transient_iterations):
        x = r * x * (1.0 - x)

    # Calculate the sum of log|f'(x)| over many iterations
    log_sum = 0.0
    for _ in range(calc_iterations):
        # Derivative of the logistic map f'(x) = r(1 - 2x)
        derivative = r * (1.0 - 2.0 * x)
        log_sum += math.log(abs(derivative))
        
        # Iterate the map to the next state
        x = r * x * (1.0 - x)

    lambda_val = log_sum / calc_iterations
    
    print(f"  - Using {calc_iterations} iterations after a {transient_iterations} transient period.")
    print(f"  - Calculated Lyapunov exponent λ ≈ {lambda_val:.7f}")
    print("-" * 50)
    
    # --- Step 3: Calculate T(n, p) ---
    print("Step 3: Calculate T(n,p) using the predictability horizon formula.")
    
    ln10 = math.log(10)
    T_np = (p * ln10) / lambda_val

    print(f"  - Formula: T(n,p) ≈ (p * ln(10)) / λ")
    print(f"  - T({n},{p}) ≈ ({p} * {ln10:.7f}) / {lambda_val:.7f}")
    print(f"  - T({n},{p}) ≈ {p * ln10:.7f} / {lambda_val:.7f}")
    print(f"  - T({n},{p}) ≈ {T_np:.7f}")
    print("-" * 50)
    
    # --- Step 4: Calculate the final answer ---
    print("Step 4: Compute the final requested value.")
    
    final_answer = math.ceil(T_np / 12)
    
    print(f"  - Formula: ceil(T(n,p) / 12)")
    print(f"  - Calculation: ceil({T_np:.7f} / 12)")
    print(f"  - Result: ceil({T_np/12:.7f})")
    print(f"  - Final Answer = {final_answer}")
    print("-" * 50)
    
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    solve_dynamical_problem()