import math

def calculate_bch_coefficients(n, k):
    """
    Calculates the number of nonzero coefficients for a given order in the BCH expansion
    using Witt's formula and prints the detailed calculation.
    """

    def get_divisors(num):
        """Finds all positive divisors of a number."""
        divs = set()
        for i in range(1, int(math.sqrt(num)) + 1):
            if num % i == 0:
                divs.add(i)
                divs.add(num // i)
        return sorted(list(divs))

    def mobius(num):
        """Calculates the Möbius function mu(n)."""
        if num == 1:
            return 1
        
        p_count = 0
        temp_n = num
        i = 2
        # Check for factors up to sqrt(temp_n)
        while i * i <= temp_n:
            if temp_n % i == 0:
                p_count += 1
                temp_n //= i
                # Check for square factors
                if temp_n % i == 0:
                    return 0
            i += 1
        
        # If temp_n is still > 1, it's a prime factor itself
        if temp_n > 1:
            p_count += 1
            
        # Return based on parity of prime factors count
        return -1 if p_count % 2 != 0 else 1

    # --- Calculation ---
    divisors_of_n = get_divisors(n)
    sum_of_terms = 0
    calculation_steps = []

    for d in divisors_of_n:
        mu_d = mobius(d)
        exponent = n // d
        term = mu_d * (k ** exponent)
        sum_of_terms += term
        calculation_steps.append({
            'd': d,
            'mu': mu_d,
            'exp': exponent,
            'term': term
        })

    result = sum_of_terms // n

    # --- Output ---
    print(f"Calculating the number of nonzero BCH coefficients of order n={n} for k={k} generators.")
    print("Using Witt's Formula: (1/n) * sum over d|n of [mu(d) * k^(n/d)]")
    
    print(f"\nStep 1: The divisors 'd' of n={n} are: {divisors_of_n}")

    print("\nStep 2: Calculate each term of the sum:")
    for step in calculation_steps:
        print(f"  For d={step['d']:<2}: term = mu({step['d']}) * {k}^({n}/{step['d']}) = ({step['mu']}) * {k}^{step['exp']} = {step['term']}")

    print(f"\nStep 3: Sum the terms and divide by n.")

    # --- Final Equation Breakdown ---
    print("The final equation is:")
    
    # Line 1: Formula with symbolic values
    line1_parts = [f"({step['mu']})*({k}^{step['exp']})" for step in calculation_steps]
    line1 = f"(1/{n}) * [ " + " + ".join(line1_parts) + " ]"
    print(line1.replace("+ (-", "- "))

    # Line 2: Evaluated terms
    line2_parts = [str(step['term']) for step in calculation_steps]
    line2 = f"= (1/{n}) * [ " + " + ".join(line2_parts) + " ]"
    print(line2.replace("+ (", "").replace(")", "").replace("+ -", "- "))

    # Line 3: Sum inside the brackets
    line3 = f"= (1/{n}) * [ {sum_of_terms} ]"
    print(line3)

    # Line 4: Final result
    line4 = f"= {result}"
    print(line4)

# Execute the calculation for the specific problem
if __name__ == "__main__":
    order_n = 10
    generators_k = 2
    calculate_bch_coefficients(order_n, generators_k)