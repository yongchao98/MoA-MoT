import numpy as np
from scipy.integrate import quad
import math

def main():
    """
    Analyzes the Borwein Integral I_n and provides data to evaluate the propositions.
    """
    print("Analysis of the Borwein Integral I_n = integral_0^inf Product_{k=1 to n} sinc(x/k) dx")
    print("A proposition P(n) states that I_n = pi/2.\n")
    print("The integral is pi/2 if and only if the condition `1 >= sum_{k=2 to n} (1/k)` holds.")
    print("-" * 80)
    print(f"{'n':<3} | {'I_n (Numerical)':<20} | {'pi/2':<20} | {'Difference (I_n-pi/2)':<25} | {'Sum Condition':<25}")
    print("-" * 80)

    # Define the sinc function as sin(x)/x
    def my_sinc(x):
        if x == 0:
            return 1.0
        return math.sin(x) / x

    # Define the integrand for I_n
    def integrand(x, n):
        prod = 1.0
        for k in range(1, n + 1):
            prod *= my_sinc(x / k)
        return prod

    current_sum = 0.0
    for n in range(1, 9):
        # Calculate the integral I_n numerically
        # The quad function is highly accurate for well-behaved integrals.
        # We integrate up to infinity.
        integral_val, error = quad(integrand, 0, np.inf, args=(n,))

        # Check the sum condition
        if n == 1:
            condition_holds = True
            sum_str = "0.0 (empty sum)"
        else:
            current_sum += 1/n
            condition_holds = (1.0 >= current_sum)
            sum_str = f"sum = {current_sum:.4f}, holds: {condition_holds}"

        # Print the results in a formatted table
        pi_half = math.pi / 2
        difference = integral_val - pi_half
        
        print(f"{n:<3} | {integral_val:<20.10f} | {pi_half:<20.10f} | {difference:<+25.15e} | {sum_str}")

if __name__ == "__main__":
    main()
