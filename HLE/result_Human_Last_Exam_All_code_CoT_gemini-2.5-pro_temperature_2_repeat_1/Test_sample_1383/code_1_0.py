import collections
from fractions import Fraction

def get_prime_factorization(n):
    """
    Returns a dictionary of the prime factors of an integer n.
    e.g., get_prime_factorization(12) -> {2: 2, 3: 1} for 12 = 2^2 * 3^1
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return dict(factors)

def format_frac(f):
    """Formats a Fraction object as a string."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"

def format_expression(coeffs_dict, name):
    """Formats a dictionary of log coefficients into a human-readable equation."""
    if not any(c != 0 for c in coeffs_dict.values()):
        return f"{name} = 0"
    
    parts = []
    # Sort keys for consistent output, e.g., ln(2), ln(3), ln(5)...
    for prime in sorted(coeffs_dict.keys()):
        coeff = coeffs_dict[prime]
        if coeff == 0:
            continue
            
        sign = "+" if coeff > 0 else "-"
        abs_coeff = abs(coeff)

        if abs_coeff == 1:
            # Hide coefficient if it's 1 or -1
            coeff_str = ""
        else:
            coeff_str = f"{format_frac(abs_coeff)} "
        
        term = f"{coeff_str}ln({prime})"
        parts.append((sign, term))
    
    # Build the final string, handling the sign of the first term
    result_str = f"{name} = "
    first_sign, first_term = parts[0]
    if first_sign == '+':
        result_str += first_term
    else:
        result_str += f"-{first_term}"

    for sign, term in parts[1:]:
        result_str += f" {sign} {term}"
        
    return result_str

def calculate_growth_coeffs(true_probs, bet_fracs, odds):
    """
    Calculates the coefficients of the log(prime) terms for a growth rate formula.
    W = sum(p_i * ln(b_i * o_i))
    """
    coeffs = collections.defaultdict(Fraction)
    for p_i, b_i, o_i in zip(true_probs, bet_fracs, odds):
        arg = b_i * o_i
        
        # log(1) is 0, so it doesn't contribute
        if arg == 1:
            continue
        
        # Decompose log(num/den) into log(num) - log(den)
        # and then into sum of log(prime)
        num_factors = get_prime_factorization(arg.numerator)
        den_factors = get_prime_factorization(arg.denominator)
        
        for prime, exponent in num_factors.items():
            coeffs[prime] += p_i * exponent
            
        for prime, exponent in den_factors.items():
            coeffs[prime] -= p_i * exponent
            
    return dict(coeffs)

def main():
    """Main function to perform the calculations and print results."""
    # Input data from the problem
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q_belief = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    odds = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]
    
    # 1. Calculate optimal growth rate W*
    # Betting fractions are the true probabilities: b_i = p_i
    coeffs_W_star = calculate_growth_coeffs(p_true, p_true, odds)
    
    # 2. Calculate achieved growth rate W with incorrect beliefs
    # Betting fractions are the incorrect probabilities: b_i = q_i
    coeffs_W = calculate_growth_coeffs(p_true, q_belief, odds)
    
    # 3. Calculate the difference Delta_W = W* - W
    all_primes = sorted(list(set(coeffs_W_star.keys()) | set(coeffs_W.keys())))
    coeffs_Delta_W = {
        p: coeffs_W_star.get(p, Fraction(0)) - coeffs_W.get(p, Fraction(0))
        for p in all_primes
    }
    
    # Print the results
    print("The optimal growth rate is:")
    print(format_expression(coeffs_W_star, "W*"))
    print("\nThe achieved growth rate with incorrect beliefs is:")
    print(format_expression(coeffs_W, "W"))
    print("\nThe decrease in the growth rate is:")
    print(format_expression(coeffs_Delta_W, "ΔW"))

if __name__ == '__main__':
    main()
