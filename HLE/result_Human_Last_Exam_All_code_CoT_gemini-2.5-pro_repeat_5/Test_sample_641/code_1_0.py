def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The number of involutions in PSU(4, q) for q=1 mod 4 is given by the formula:
    # (q^8 - q^7 + 2*q^6 + q^4 + q^2) / 2
    # We will calculate each term and then the final result.

    q2 = q**2
    q4 = q**4
    q6 = q**6
    q7 = q**7
    q8 = q**8

    numerator = q8 - q7 + 2 * q6 + q4 + q2
    result = numerator // 2

    print("The number of involutions in PSU(4, 997) is calculated using the formula:")
    print("(q^8 - q^7 + 2*q^6 + q^4 + q^2) / 2, where q = 997.")
    print("\nStep 1: Calculate the powers of q:")
    print(f"q^2 = {q}^2 = {q2}")
    print(f"q^4 = {q}^4 = {q4}")
    print(f"q^6 = {q}^6 = {q6}")
    print(f"q^7 = {q}^7 = {q7}")
    print(f"q^8 = {q}^8 = {q8}")

    print("\nStep 2: Plug these values into the numerator of the formula:")
    print(f"Numerator = q^8 - q^7 + 2*q^6 + q^4 + q^2")
    print(f"Numerator = {q8} - {q7} + 2*{q6} + {q4} + {q2}")
    print(f"Numerator = {numerator}")

    print("\nStep 3: Divide by 2 to get the final answer:")
    print(f"Number of involutions = {numerator} / 2")
    print(f"Number of involutions = {result}")

solve()