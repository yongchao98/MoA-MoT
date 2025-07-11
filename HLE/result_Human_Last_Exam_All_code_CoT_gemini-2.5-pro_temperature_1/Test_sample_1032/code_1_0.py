def solve_horse_area():
    """
    Calculates the reachable area for the horse.
    """
    # 1. Define the parameters
    # The rope has a taxi-cab length of L = 7/2
    numerator = 7
    denominator = 2
    rope_length_L = numerator / denominator

    # The house has an area of 3 unit squares
    house_area = 3

    # 2. Calculate the area of the reachable diamond without the house
    # The area of a taxi-cab circle (diamond) with radius L is 2 * L^2
    diamond_area = 2 * (rope_length_L ** 2)

    # 3. Calculate the final reachable area
    # The house is located at a corner of the origin, within the diamond.
    # The rope cannot enter the house.
    # Because the rope is tied to a vertex of the house and the house does not
    # create any "shadows" for straight-line paths from the origin,
    # the reachable area is the diamond's area minus the house's area.
    reachable_area = diamond_area - house_area

    # 4. Print the final equation with all its numbers
    # The required numbers are: 2, 7, 2, 2, 3, and the result 21.5
    # The equation is: 2 * (7/2)^2 - 3 = 21.5
    val_1 = 2
    val_2 = numerator
    val_3 = denominator
    val_4 = 2  # The exponent
    val_5 = house_area
    result = reachable_area
    
    print("The final area is the area of the taxi-cab diamond minus the area of the house.")
    print("The calculation is:")
    print(f"{val_1} * ({val_2}/{val_3})^{val_4} - {val_5} = {result}")

solve_horse_area()