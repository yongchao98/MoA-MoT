import math

def solve_mixing_problem():
    """
    Calculates the number of repetitions needed for the liquids in two cups to become "the same".
    """
    # Fix parameters epsilon and delta as per the problem description.
    # You can change these values to see the result for different parameters.
    epsilon = 0.01
    delta = 0.1

    print("This script calculates the number of repetitions (N) for the mixing process.")
    print(f"Given parameters: ε = {epsilon}, δ = {delta}\n")

    # --- Step-by-step explanation ---
    print("--- Derivation of the Formula ---")
    print("Step 1: Model the mixing process.")
    print("Let r_A(n) be the concentration of red liquid in cup A after n repetitions.")
    print("Assuming a simultaneous swap of a δ-fraction of liquid, the total volume in each cup remains 1.")
    print("The recurrence relation is: r_A(n+1) = r_A(n) - δ*r_A(n) + δ*r_B(n).")
    print("Using r_A(n) + r_B(n) = 1, we get: r_A(n+1) = (1-2δ)*r_A(n) + δ.\n")

    print("Step 2: Solve the recurrence relation.")
    print("With the initial condition r_A(0)=1, the solution is: r_A(n) = 0.5 + 0.5 * (1-2δ)^n.\n")

    print("Step 3: Define the stopping condition.")
    print("The liquids are 'the same' if |r_A(n) - r_B(n)| <= ε.")
    print("The difference is |r_A(n) - (1 - r_A(n))| = |2*r_A(n) - 1| = |(1-2δ)^n|.")
    print(f"So we need the smallest integer n such that |(1 - 2*δ)|^n <= ε.\n")

    print("Step 4: Solve for n.")
    print("Taking logarithms, we get: n * log(|1 - 2δ|) <= log(ε).")
    print("Since |1 - 2δ| < 1 (for δ in (0,1) and δ != 0.5), its log is negative, so we flip the inequality:")
    print("n >= log(ε) / log(|1 - 2δ|).")
    print("N is the smallest integer n, which is the ceiling of this expression.\n")

    # --- Calculation ---
    print("--- Calculation for Given Values ---")
    
    # Handle edge cases
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon (ε) and delta (δ) must be between 0 and 1.")
        return

    if delta == 0.5:
        N = 1
        print("Special case: For δ = 0.5, perfect mixing occurs in exactly one step.")
        print(f"\nFinal Answer: N = {N}")
    else:
        log_epsilon = math.log(epsilon)
        term = 1 - 2 * delta
        abs_term = abs(term)
        log_abs_term = math.log(abs_term)
        
        value = log_epsilon / log_abs_term
        # The number of repetitions must be a non-negative integer.
        # max(0,...) handles cases where epsilon >= 1, resulting in a non-positive N.
        N = math.ceil(max(0, value))

        print("Final Equation:")
        print(f"N = ceil(log({epsilon}) / log(|1 - 2*{delta}|))")
        print(f"N = ceil({log_epsilon:.4f} / log(|{term:.2f}|))")
        print(f"N = ceil({log_epsilon:.4f} / log({abs_term:.2f}))")
        print(f"N = ceil({log_epsilon:.4f} / {log_abs_term:.4f})")
        print(f"N = ceil({value:.4f})")
        print(f"\nFinal Answer: N = {N}")

if __name__ == '__main__':
    solve_mixing_problem()