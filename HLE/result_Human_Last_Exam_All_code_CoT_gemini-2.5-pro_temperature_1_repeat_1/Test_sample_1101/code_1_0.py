import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(3/2+n) / Gamma(3/2+n-m)
    """
    if n < 0:
        return 0
    
    total_sum = 0
    
    # Pre-calculate the constant part for a given n
    try:
        gamma_n_plus_3_2 = math.gamma(1.5 + n)
    except ValueError:
        print(f"Error: Gamma function not defined for {1.5 + n}")
        return None

    print(f"Calculating S_{n} with n = {n}:")
    print("Equation: S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n) / Gamma(1.5+n-m)")
    print("--------------------------------------------------")
    
    for m in range(n + 1):
        # Binomial coefficient C(n, m)
        comb = math.comb(n, m)
        
        # Gamma term
        try:
            gamma_term_val = math.gamma(1.5 + n - m)
            if gamma_term_val == 0:
                # This case is unlikely for the given arguments but good practice
                term_val = 0
            else:
                gamma_ratio = gamma_n_plus_3_2 / gamma_term_val
        except ValueError:
            # This can happen if 1.5 + n - m is a non-positive integer
            gamma_ratio = 0

        term = ((-1)**m) * comb * gamma_ratio
        
        # Printing each term of the final equation as requested
        sign = "-" if m % 2 != 0 else "+"
        if m == 0:
          sign = " "
        print(f"m={m}: {sign} {comb} * {gamma_ratio:.4f} = {term:.4f}")
        
        total_sum += term
        
    return total_sum

def main():
    n_val = 10
    
    s_n = calculate_sum(n_val)
    
    if s_n is not None:
      print("--------------------------------------------------")
      print(f"Final Sum S_{n_val} = {s_n}")
      
      # Asymptotic behavior check
      if n_val > 0:
          f_n = math.exp(2 * math.sqrt(n_val))
          ratio = s_n / f_n
          print(f"The proposed bounding function is f(n) = e^(2*sqrt(n))")
          print(f"f({n_val}) = {f_n}")
          print(f"The ratio S_{n_val} / f({n_val}) = {ratio}")
          
          # The theoretical constant for the ratio
          const = 1 / (2 * math.sqrt(math.pi) * math.sqrt(math.e))
          print(f"Asymptotic theory suggests this ratio should approach C = {const}")

if __name__ == "__main__":
    main()
