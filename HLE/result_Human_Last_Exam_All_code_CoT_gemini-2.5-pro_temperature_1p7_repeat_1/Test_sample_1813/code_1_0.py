import math

def get_continued_fraction_coeffs(p, q):
    """
    Computes the continued fraction coefficients [a_1, a_2, ..., a_n] for a rational number p/q,
    where it is assumed 0 < p/q < 1. The first coefficient a_0=0 is implicit.
    """
    if p == 0:
        return []
    coeffs = []
    # Make copies to not modify the original numbers
    num, den = p, q
    while num > 0:
        a = den // num
        coeffs.append(a)
        den, num = num, den % num
    return coeffs

def compute_generalized_markov_number(p, q, coeffs):
    """
    Computes the generalized Markov number using the continued fraction coefficients.
    The recurrence relation used is x_{i+1} = a_i * x_i - x_{i-1} with x_0=0, x_1=1.
    """
    print(f"Step 2: Compute the generalized Markov number m_{p}/{q}")
    
    # Initial conditions for the recurrence
    x_prev = 0  # Represents x_0
    x_curr = 1  # Represents x_1
    
    print("Applying the recurrence relation: x_{i+1} = a_i * x_i - x_{i-1}")
    print(f"Initial values: x_0 = {x_prev}, x_1 = {x_curr}")

    # Loop through the coefficients to calculate the terms of the sequence
    for i, a in enumerate(coeffs):
        x_next = a * x_curr - x_prev
        print(f"For a_{i+1}={a}: x_{i+2} = {a} * x_{i+1} - x_i = {a} * {x_curr} - {x_prev} = {x_next}")
        x_prev, x_curr = x_curr, x_next
        
    return x_curr

def get_continued_fraction_string(n):
    """
    Computes the standard continued fraction representation for an integer n.
    For an integer, the representation is simply [n].
    """
    return f"[{n}]"

def main():
    """
    Main function to orchestrate the computation and printing of the result.
    """
    p, q = 4, 7

    print(f"Task: Compute the continued fraction associated with the generalized Markov number m_{p}/{q}.\n")
    
    # Step 1: Compute continued fraction of p/q
    print(f"Step 1: Compute the continued fraction of {p}/{q}")
    cf_coeffs = get_continued_fraction_coeffs(p, q)
    # The standard notation for p/q < 1 is [0; a_1, a_2, ...]
    cf_str_pq = f"[0; {', '.join(map(str, cf_coeffs))}]"
    print(f"The continued fraction of {p}/{q} is {cf_str_pq}.")
    print("-" * 30)

    # Step 2: Compute the generalized Markov number
    markov_number = compute_generalized_markov_number(p, q, cf_coeffs)
    print(f"\nThe generalized Markov number m_{p}/{q} is: {markov_number}")
    print("-" * 30)

    # Step 3: Compute the continued fraction of the resulting Markov number
    print(f"Step 3: Compute the continued fraction of the result ({markov_number})")
    final_cf_str = get_continued_fraction_string(markov_number)
    print(f"The continued fraction of {markov_number} is {final_cf_str}.")
    print("-" * 30)
    
    # Final Answer
    print("The final result presents the generalized Markov number and its continued fraction representation in an equation.")
    print(f"Final Equation: {markov_number} = {final_cf_str}")

if __name__ == "__main__":
    main()
