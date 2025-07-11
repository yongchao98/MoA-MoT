def solve_pearl_riddle():
    """
    This function solves the pearl necklace riddle step-by-step.
    """
    # Step 1: Calculate the number of pearls remaining on the string.
    # "seven shy of eleven times eleven"
    eleven = 11
    seven = 7
    remaining_on_string = (eleven * eleven) - seven

    print(f"First, let's find the number of pearls that remained on the string.")
    print(f"This is 'seven shy of eleven times eleven', which calculates to: ({eleven} * {eleven}) - {seven} = {remaining_on_string}\n")

    # Step 2: Solve for the total number of pearls ('x').
    # The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remaining_on_string
    # This can be rearranged to: x - (1/6 + 1/5 + 1/3 + 1/10)x = remaining_on_string
    # The sum of the fractions is 24/30, which simplifies to 4/5.
    # So, x * (1 - 4/5) = remaining_on_string  =>  x * (1/5) = remaining_on_string
    # Therefore, x = remaining_on_string * 5
    total_pearls = remaining_on_string * 5

    print("Next, we solve for the total number of pearls on the necklace.")
    print("The full equation representing the state of the pearls is:")
    
    # To show each number in the final equation, we calculate each part
    floor_pearls = total_pearls * (1/6)
    bed_pearls = total_pearls * (1/5)
    woman_pearls = total_pearls * (1/3)
    lover_pearls = total_pearls * (1/10)
    
    print(f"Total Pearls ({int(total_pearls)}) = On Floor ({int(floor_pearls)}) + On Bed ({int(bed_pearls)}) + With Woman ({int(woman_pearls)}) + With Lover ({int(lover_pearls)}) + On String ({int(remaining_on_string)})")
    print(f"So, there were {int(total_pearls)} pearls altogether.\n")

    # Step 3: Calculate how many more pearls are needed.
    fallen_pearls = total_pearls - remaining_on_string
    found_pearls = fallen_pearls / 3
    needed_pearls = fallen_pearls - found_pearls

    print("Finally, let's calculate how many more pearls they need.")
    print(f"The number of pearls that fell off is the total minus those remaining on the string: {int(total_pearls)} - {int(remaining_on_string)} = {int(fallen_pearls)}.")
    print(f"They manage to find 1/3rd of the fallen pearls: {int(fallen_pearls)} / 3 = {int(found_pearls)}.")
    print(f"The number of pearls they still need is the number that fell minus the number they found: {int(fallen_pearls)} - {int(found_pearls)} = {int(needed_pearls)}.")
    print(f"\nTherefore, they are going to need {int(needed_pearls)} more pearls.")

solve_pearl_riddle()