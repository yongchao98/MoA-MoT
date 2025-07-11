def number_of_feynman_diagrams(n):
    """
    Calculates a(n), the number of non-vanishing Feynman diagrams of order 2n.
    This function uses a recurrence relation valid for small n.
    """
    if n == 0:
        return 1
    if n == 1:
        return 1
    
    # Calculate a(n-1) and a(n-2) recursively
    val_n_minus_1 = number_of_feynman_diagrams(n-1)
    val_n_minus_2 = number_of_feynman_diagrams(n-2)

    # Apply the recurrence relation
    # a(n) = (5n - 7) * a(n-1) + 4 * (n-1)^2 * a(n-2)
    coeff1 = 5 * n - 7
    coeff2 = 4 * (n - 1)**2
    
    result = coeff1 * val_n_minus_1 + coeff2 * val_n_minus_2
    
    # For the final answer, print the equation details
    if n == 3:
        print(f"Calculating a(3) using the recurrence a(n) = (5n - 7)*a(n-1) + 4*(n-1)^2*a(n-2):")
        a2 = number_of_feynman_diagrams(2)
        a1 = number_of_feynman_diagrams(1)
        print(f"a(3) = (5*3 - 7)*a(2) + 4*(3-1)^2*a(1)")
        print(f"a(3) = ({coeff1})*{a2} + ({coeff2})*{a1}")
        print(f"a(3) = {coeff1 * a2} + {coeff2 * a1}")
        print(f"a(3) = {result}")

    return result

# Calculate and print a(3)
final_answer = number_of_feynman_diagrams(3)
print(f"\nThe value of a(3) is: {final_answer}")