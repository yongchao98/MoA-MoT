import math

def print_rates(T):
    """Prints the values of different convergence rate functions for a given T."""
    # Option A: 1/T
    rate_a = 1 / T
    
    # Option B: 1/T^(2/3)
    rate_b = 1 / (T**(2/3))
    
    # Option C: 1/T^(1/2)
    rate_c = 1 / math.sqrt(T)
    
    # Derived Rate (for option D): log(T)/sqrt(T)
    rate_d = math.log(T) / math.sqrt(T)
    
    print(f"For T = {T:.0e}:")
    print(f"  Rate A (1/T)          = {rate_a:.2e}")
    print(f"  Rate B (1/T^(2/3))      = {rate_b:.2e}")
    print(f"  Rate C (1/sqrt(T))      = {rate_c:.2e}")
    print(f"  Rate D (log(T)/sqrt(T)) = {rate_d:.2e}")
    print("-" * 20)

if __name__ == "__main__":
    for power in [3, 6, 9, 12]:
        T_val = 10**power
        print_rates(T_val)
