def solve_minecraft_riddle():
    """
    Solves the Minecraft riddle by identifying the first obtainable block in the given scenario.

    The scenario describes an impossible starting condition in vanilla Minecraft unless
    a "Bonus Chest" is enabled at world creation. This chest provides the initial
    resources needed to advance.

    The most fundamental block for advancement found in a bonus chest is the Wood Log,
    as it unlocks the entire crafting tree.
    """

    # The first block essential for game progression.
    first_block = "Wood Log"

    # The crafting recipe for a Wood Log can be seen as an equation.
    # The numbers involved are 1 (for the log) and 4 (for the planks).
    input_block_count = 1
    output_plank_count = 4

    print(f"The first block that must be obtained to advance in the game is the: {first_block}")
    print("\nThis block is the key to all further progress.")
    print("It enables the game's first and most important crafting equation:")
    print(f"{input_block_count} {first_block} -> {output_plank_count} Wood Planks")

solve_minecraft_riddle()