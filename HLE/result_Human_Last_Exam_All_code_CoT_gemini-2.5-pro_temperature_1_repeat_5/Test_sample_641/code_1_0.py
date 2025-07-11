def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    The formula for the number of involutions in PSU(4, q) for q = 1 (mod 4)
    is N = q^4 * (q^4 + q^2 + 1).
    """
    q = 997

    # Using Python's arbitrary-precision integers for large number calculations
    q2 = q * q
    q4 = q2 * q2

    # The final formula
    num_involutions = q4 * (q4 + q2 + 1)

    # Print the explanation and the final equation with numbers
    print("The number of involutions in PSU(4, 997) is calculated using the formula:")
    print("N = q^4 * (q^4 + q^2 + 1), for q = 997\n")
    print("Step-by-step calculation:")
    print(f"Let q = {q}")
    print(f"q^2 = {q}^2 = {q2}")
    print(f"q^4 = ({q2})^2 = {q4}")
    
    sum_term = q4 + q2 + 1
    
    print("\nSubstitute these values into the formula:")
    print(f"N = {q4} * ({q4} + {q2} + 1)")
    print(f"N = {q4} * {sum_term}")
    print(f"N = {num_involutions}")

solve()