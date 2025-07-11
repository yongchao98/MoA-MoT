def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The formula for the number of involutions is (q**4 * (q**2 + 1) * (q**2 - q + 1)) / 2
    term1 = q**4
    term2 = q**2 + 1
    term3 = q**2 - q + 1
    
    result = (term1 * term2 * term3) // 2

    # Print the explanation and the step-by-step calculation
    print("The number of involutions in PSU(4, 997) is given by the formula:")
    print("(q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2, where q = 997.\n")
    print("Step-by-step calculation:")
    print(f"q = {q}")
    print(f"q^2 = {q**2}")
    print(f"q^4 = {term1}")
    print(f"q^2 + 1 = {term2}")
    print(f"q^2 - q + 1 = {term3}\n")
    
    print("The full equation is:")
    print(f"({term1} * {term2} * {term3}) / 2\n")

    print("Final result:")
    print(result)

solve()