def solve_involutions():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    # Define the prime field characteristic q
    q = 997

    # The number of involutions in PSU(4, q) for q = 1 (mod 4) is given by the formula:
    # N = (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2
    # We use integer arithmetic throughout.

    # Calculate the intermediate terms
    q_2 = q * q
    q_4 = q_2 * q_2
    term1 = q_2 - q + 1
    term2 = q_2 + 1

    # The numerator of the formula
    numerator = q_4 * term1 * term2

    # The final result, using integer division
    num_involutions = numerator // 2

    # Print the explanation and the numbers involved in the calculation
    print(f"To find the number of involutions in PSU(4, {q}), we use the formula:")
    print("N = (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2\n")
    print(f"First, we calculate the components with q = {q}:")
    print(f"q^4 = {q}**4 = {q_4}")
    print(f"q^2 - q + 1 = {q}**2 - {q} + 1 = {q_2} - {q} + 1 = {term1}")
    print(f"q^2 + 1 = {q}**2 + 1 = {q_2} + 1 = {term2}\n")

    print("The final equation with these numbers is:")
    print(f"N = ({q_4} * {term1} * {term2}) / 2\n")

    print("The final result is:")
    print(num_involutions)

solve_involutions()