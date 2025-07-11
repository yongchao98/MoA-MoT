import io
import sys

# Set stdout to use UTF-8 encoding
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def solve_riddle():
    """
    This function explains the solution to the Chinese character riddle.
    """
    # The components of the riddle
    wooden_stick = "木 (mù)"
    box_and_ladder = "覀 (yà)"
    
    # The final character
    final_character = "栗 (lì)"

    print("Let's solve the riddle step-by-step by deconstructing the character:")
    print("-" * 60)
    
    # Explain the first part of the riddle
    print("1. 'A wooden stick' clearly points to the character radical for 'wood':")
    print(f"   - Component 1: {wooden_stick}")
    
    # Explain the second part of the riddle
    print("\n2. 'hanging a square box, a ladder placed in the center' describes the top part.")
    print("   This is the component '覀', which can be visualized as:")
    print("   - A 'square box' (the 冂 shape)")
    print("   - With a 'ladder' (the two vertical strokes '丨丨') inside it.")
    print(f"   - Component 2: {box_and_ladder}")
    
    # Combine the parts
    print("\n3. Combining these components, with '覀' on top of '木', gives us the answer.")
    print(f"   {box_and_ladder} + {wooden_stick} => {final_character}")
    
    print("-" * 60)
    print(f"The character is: {final_character.split()[0]}")

solve_riddle()