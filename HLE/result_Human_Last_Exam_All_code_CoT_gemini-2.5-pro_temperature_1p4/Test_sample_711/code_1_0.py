def solve_minecraft_riddle():
    """
    Solves the riddle of the first obtainable block in a constrained Minecraft world.
    """
    # Define the components of the solution based on the riddle's logic.
    mob_source = "Witch"
    required_item = "Glowstone Dust"
    required_item_count = 4
    resulting_block = "Glowstone Block"
    crafting_grid = "2x2 inventory grid"

    # Explain the deduction process.
    print("In this scenario, you cannot mine or find blocks. The only way to get new items is from mob drops.")
    print(f"A {mob_source} is the only mob that drops an item which can be crafted into a block in the {crafting_grid}.")
    print(f"By defeating Witches, you can collect {required_item}.")
    
    # Print the final crafting equation, showing each number and item name.
    print("\nThe first obtainable block is crafted as follows:")
    print(f"{required_item_count} {required_item} -> 1 {resulting_block}")

solve_minecraft_riddle()