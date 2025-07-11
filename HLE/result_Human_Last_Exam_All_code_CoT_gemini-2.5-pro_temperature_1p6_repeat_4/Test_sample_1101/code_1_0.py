import math

def calculate_sum(n_max):
    """
    Calculates the sum S_n for n from 0 to n_max.
    S_n = sum_{m=0 to n} (-1)^m * C(n, m) * Product_{j=1 to m} (n + 3/2 - j)
    """
    
    print("n\tS_n\t\t|S_n| / (4^n)")
    print("-" * 40)
    
    for n in range(n_max + 1):
        total_sum = 0.0
        for m in range(n + 1):
            # Calculate binomial coefficient C(n, m)
            try:
                n_choose_m = math.comb(n, m)
            except ValueError:
                n_choose_m = 0

            # Calculate the product term (n + 1/2)_m
            pochhammer = 1.0
            for i in range(m):
                pochhammer *= (n + 0.5 - i)
            
            term = ((-1)**m) * n_choose_m * pochhammer
            total_sum += term
        
        # Analyze the asymptotic behavior
        if n == 0:
            ratio = abs(total_sum)
        else:
            ratio = abs(total_sum) / (4**n)

        print(f"{n}\t{total_sum:10.4f}\t{ratio:10.6f}")

# Calculate for n up to 15 to observe the trend
calculate_sum(15)
