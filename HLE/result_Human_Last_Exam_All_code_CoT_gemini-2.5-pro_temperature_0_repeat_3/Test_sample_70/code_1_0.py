def solve_flag_rank():
    """
    Identifies African nations whose flags have the same linear algebraic rank as the flag of Denmark.

    The rank of the Danish flag is 2. We look for African flags that are also rank 2.
    This occurs when a flag's design is composed of exactly two linearly independent row or column patterns.
    """

    # The rank of the Danish flag is 2. The final equation for any matching flag is rank(Flag) = 2.
    # The number in this equation is 2.
    denmark_rank = 2

    print(f"The rank of the flag of Denmark is {denmark_rank}.")
    print("The final equation for a flag to match is: rank(Flag) = 2\n")

    # List of African countries whose flags have a rank of 2
    rank_2_african_flags = ["Benin", "Madagascar", "Nigeria"]

    print("The flags of the following African nations have the same rank:")
    for country in rank_2_african_flags:
        # For each country, we output the final equation and its components.
        # The equation is rank(Country) = 2. The number is 2.
        print(f"rank({country}) = {denmark_rank}")

    print("\nFinal list of countries:")
    for country in rank_2_african_flags:
        print(country)

solve_flag_rank()