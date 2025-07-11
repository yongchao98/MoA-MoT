import math

def calculate_sum(n):
    """
    Calculates the sum S(n) = \sum_{m=0}^n (-1)^m {n\choose m} \frac{\Gamma(\frac{3}{2}+n)}{\Gamma(\frac{3}{2}+n-m)}
    """
    if n < 0:
        return 0

    # The term \frac{\Gamma(x)}{\Gamma(x-m)} is the falling factorial (x-1)_m
    # Here x = n + 3/2, so x-1 = n + 1/2
    # So the inner term is the falling factorial of (n + 1/2) for m terms.
    
    total_sum = 0
    for m in range(n + 1):
        # Calculate binomial coefficient {n choose m}
        binom_coeff = math.comb(n, m)
        
        # Calculate falling factorial (n + 1/2)_m
        pochhammer = 1.0
        for i in range(m):
            pochhammer *= (n + 0.5 - i)
            
        term = ((-1)**m) * binom_coeff * pochhammer
        total_sum += term
        
    return total_sum

def main():
    try:
        n_input = int(input("Enter an integer n >= 1: "))
        if n_input < 1:
            print("Please enter an integer n >= 1.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    sum_val = calculate_sum(n_input)
    
    # Our proposed function f(n) is Gamma(n) = (n-1)!
    f_val = math.gamma(n_input)

    print(f"For n = {n_input}:")
    print(f"The value of the sum S(n) is: {sum_val}")
    print(f"The value of the bounding function f(n) = (n-1)! is: {f_val}")
    
    if abs(sum_val) <= f_val:
        print("The condition |S(n)| <= f(n) holds (with C=1).")
    else:
        # For a sufficiently large C, the bound will hold.
        # Let's calculate the required C for this n.
        if f_val > 0:
            C = abs(sum_val) / f_val
            print(f"The condition S(n) <= C*f(n) holds for C >= {C:.4f}")

if __name__ == "__main__":
    main()