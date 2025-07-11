# Given an integer k >= 2.
# You can change the value of k to compute the limit for a different k.
k = 4

if not isinstance(k, int) or k < 2:
    print("Error: k must be an integer greater than or equal to 2.")
else:
    # The derived formula for the limit is k / (2*k - 1)
    numerator = k
    denominator = 2 * k - 1
    result = numerator / denominator

    # Output the final equation and the result, showing each number.
    print(f"The formula for the limit is: k / (2*k - 1)")
    print(f"For the given k = {k}, the calculation is:")
    print(f"{k} / (2 * {k} - 1) = {numerator} / {denominator} = {result}")
