import math

def get_cf_quadratic(p0, q0, d):
    """
    Computes the continued fraction of a quadratic irrational of the form (p0 + sqrt(d)) / q0.
    
    Args:
        p0 (int): Initial numerator part.
        q0 (int): Initial denominator.
        d (int): The number under the square root.

    Returns:
        A list of coefficients for the continued fraction.
    """
    
    # Check that d is not a perfect square
    sqrt_d = math.sqrt(d)
    if sqrt_d == int(sqrt_d):
        raise ValueError("d must not be a perfect square.")

    coefficients = []
    seen_states = {}
    
    p = p0
    q = q0
    
    # Limit to prevent infinite loops in case of error, though it should terminate for quadratic irrationals.
    for _ in range(100):
        state = (p, q)
        if state in seen_states:
            # Period detected
            start_index = seen_states[state]
            non_periodic_part = coefficients[:start_index]
            periodic_part = coefficients[start_index:]
            
            print(f"The continued fraction is periodic.")
            # Format the output for the final equation form [a0; a1, ..., (a_k, ...)]
            if non_periodic_part:
                print(f"Non-periodic part: {non_periodic_part}")
            print(f"Periodic part: {periodic_part}")
            
            # Per instruction: "output each number in the final equation!"
            print("\nThe numbers in the final continued fraction equation are:")
            all_coeffs = non_periodic_part + periodic_part
            for coeff in all_coeffs:
                print(coeff)
            return

        seen_states[state] = len(coefficients)
        
        # Calculate the next coefficient
        # a_i = floor( (p_i + sqrt(d)) / q_i )
        a = int((p + sqrt_d) / q)
        coefficients.append(a)
        
        # Update p and q for the next iteration
        # p_{i+1} = a_i * q_i - p_i
        p_next = a * q - p
        # q_{i+1} = (d - p_{i+1}^2) / q_i
        q_next = (d - p_next**2) // q
        
        p = p_next
        q = q_next

    print("Could not find period within 100 iterations.")


if __name__ == "__main__":
    # The number is lambda = 4 + 2*sqrt(2) = 4 + sqrt(8)
    # We represent it as (p0 + sqrt(d)) / q0
    p0 = 4
    d = 8  # 2*sqrt(2) = sqrt(8)
    q0 = 1
    
    print(f"Computing the continued fraction associated with m_4/7.")
    print(f"This is the continued fraction for lambda_4/7 = 4 + 2*sqrt(2) = (4 + sqrt(8))/1.")
    print("-" * 20)
    
    get_cf_quadratic(p0, q0, d)
