def solve_minecraft_riddle():
    """
    This function determines the first block obtainable in the described
    Minecraft scenario and prints the logic and crafting equation.
    """
    
    # In a void world on an obsidian platform, the only renewable resources
    # are from hostile mob drops (Zombies, Skeletons, Creepers).
    
    # We need to craft a block. The two main possibilities from mob drops are:
    # 1. Iron Block: from Iron Ingots (rare drop from Zombies)
    # 2. Bone Block: from Bones (common drop from Skeletons)
    
    # The Bone Block is far more reliable to obtain first.
    
    item_name = "Bone"
    block_name = "Bone Block"
    
    # Crafting recipe numbers
    items_needed_for_craft = 9
    blocks_produced_from_craft = 1
    
    print("The first obtainable block is the Bone Block.")
    print("This is because Skeletons, which commonly spawn at night, drop bones.")
    print("These bones can be crafted into a block in the inventory.")
    print("\nThe crafting equation is as follows:")
    
    # Final equation outputting each number
    print(f"{items_needed_for_craft} {item_name}s -> {blocks_produced_from_craft} {block_name}")

solve_minecraft_riddle()