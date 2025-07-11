def solve_vector_field_zeros():
    """
    Calculates the minimum number of zeros for a vector field on a compact
    manifold M with a non-empty boundary, based on the Euler characteristics
    of M and its boundary.
    """
    # Example: A 2-dimensional disk (M = D^2)
    # The Euler characteristic of the disk is 1.
    chi_M = 1
    # The boundary of the disk is a circle (∂M = S^1).
    # The Euler characteristic of a circle is 0.
    chi_partial_M = 0

    print(f"Consider a manifold M with Euler characteristic chi(M) = {chi_M}.")
    print(f"Its boundary ∂M has Euler characteristic chi(∂M) = {chi_partial_M}.")
    print("\nThe least number of zeros of a vector field on M is 0 if 2*chi(M) - chi(∂M) = 0, and 1 otherwise.")
    print("\nCalculating the test value:")

    # Calculate the value of the expression 2*χ(M) - χ(∂M)
    test_value = 2 * chi_M - chi_partial_M

    # Output the equation with the numbers plugged in
    print(f"2 * {chi_M} - {chi_partial_M} = {test_value}")

    # Determine the result based on whether the test value is zero
    if test_value == 0:
        result = 0
        print("\nSince the result of the expression is 0, M admits a nowhere-zero vector field.")
    else:
        result = 1
        print("\nSince the result of the expression is not 0, any vector field on M must have at least one zero.")

    print(f"The least number of zeros is: {result}")
    
    # Return the final integer answer for the example case.
    return result

# Execute the function and print the final answer
final_answer = solve_vector_field_zeros()
# The final answer is wrapped according to the format requirement.
# For the example of a 2D-disk, the answer is 1.
print(f"\n<<< {final_answer} >>>")
