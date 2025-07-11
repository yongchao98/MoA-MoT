def solve_probability():
    """
    Calculates the probability P_m for a given m.
    The user can change the value of m inside the function.
    """
    # The user can change the value of m here
    m = 10

    # According to our analysis, the number of favorable pairs (i, j)
    # where the sequence is (i,j)-divisible is m + 2.
    num_favorable = m + 2

    # The total number of ways to choose two distinct indices i and j from
    # a set of 4m+2 is (4m+2) choose 2.
    # C(4m+2, 2) = (4m+2)(4m+1)/2 = (2m+1)(4m+1)
    term1_total = 2 * m + 1
    term2_total = 4 * m + 1
    num_total = term1_total * term2_total

    # The probability P_m is the ratio of favorable outcomes to the total number of outcomes.
    print(f"Given m = {m}")
    print("The formula for the probability is P_m = (m + 2) / ((2*m + 1) * (4*m + 1))")
    print("--------------------------------------------------")
    print("Calculating the result:")
    
    # We print each number in the final equation as requested
    print(f"P_{m} = ({m} + 2) / ((2*{m} + 1) * (4*{m} + 1))")
    print(f"P_{m} = {num_favorable} / ({term1_total} * {term2_total})")
    print(f"P_{m} = {num_favorable} / {num_total}")

    # Also printing the result as a floating-point number
    probability_decimal = num_favorable / num_total
    print(f"The decimal value of the probability is approximately: {probability_decimal}")

solve_probability()
