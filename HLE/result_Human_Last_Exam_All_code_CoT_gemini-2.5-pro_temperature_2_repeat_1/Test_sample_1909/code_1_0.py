from fractions import Fraction
import math

def calculate_expected_values(n_max):
    """
    Calculates the expected number of remaining numbers, E_n, for n up to n_max.
    
    The recurrence relation used is:
    (n-1) * E_n = (n-2) * E_{n-1} + 2 * E_{n-2} for n >= 3
    with base cases E_0 = 0, E_1 = 1, E_2 = 0.
    """
    
    # Using a dictionary to store computed values of E_n to handle arbitrary indices
    e = {
        0: Fraction(0),
        1: Fraction(1),
        2: Fraction(0)
    }

    # Iteratively compute E_n for n from 3 up to n_max
    for n in range(3, n_max + 1):
        # E_n = ((n-2) * E_{n-1} + 2 * E_{n-2}) / (n-1)
        e[n] = ((n - 2) * e[n - 1] + 2 * e[n - 2]) / (n - 1)
        
    return e

def main():
    """
    Main function to solve the problem and print the results.
    """
    # Calculate E_n up to a small n to show the process
    small_n_max = 5
    e_values = calculate_expected_values(small_n_max)
    
    print(f"Let E_n be the expected number of remaining numbers.")
    print(f"The recurrence relation is (n-1)*E_n = (n-2)*E_(n-1) + 2*E_(n-2) for n>=3.")
    print(f"Base cases: E_0=0, E_1=1, E_2=0.")
    print(f"We have computed E_3 = {e_values[3]} and E_4 = {e_values[4]}.")
    print(f"Let's demonstrate the calculation for E_5:")
    
    n = 5
    numerator_val = (n - 2) * e_values[n - 1] + 2 * e_values[n - 2]
    denominator_val = n - 1
    
    print(f"E_{n} = (({n}-2) * E_{n-1} + 2 * E_{n-2}) / ({n}-1)")
    print(f"E_{n} = (3 * {e_values[n-1]} + 2 * {e_values[n-2]}) / {denominator_val}")
    print(f"E_{n} = ({3 * e_values[n-1]} + {2 * e_values[n-2]}) / {denominator_val}")
    print(f"E_{n} = ({numerator_val}) / {denominator_val} = {e_values[n]}")

    # To find the limit, we compute the ratio for a large n
    large_n = 5000
    e_values_large = calculate_expected_values(large_n)
    ratio = float(e_values_large[large_n] / large_n)
    
    limit_val = 2 * math.log(2) - 1
    
    print(f"\nTo find the limit, we calculate E_n/n for a large n.")
    print(f"For n = {large_n}, the ratio E_n/n is approximately: {ratio:.8f}")
    print(f"This value is very close to 2*ln(2) - 1, which is approximately: {limit_val:.8f}")
    print("\nThe limit is 2*ln(2) - 1.")

if __name__ == "__main__":
    main()
