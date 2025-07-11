import sys

def solve_riddle():
    """
    This function explains the logical steps to solve the riddle
    about the English poet's surname found in a Russian text about Vienna.
    """

    print("Step 1: The riddle points to the 'wide boulevards of Vienna'. Vienna's boulevards, part of the Ringstrasse, are famous for being lush and green.")

    print("\nStep 2: In the Russian language, a very common way to describe something green, like a leafy boulevard, is with the adjective 'зелёный' (zelyonyy).")
    
    description_in_russian = "зелёный бульвар"
    print(f"   The full phrase for 'green boulevard' would be: '{description_in_russian}'.")

    print("\nStep 3: The key to the riddle is a phonetic pun. The Russian word 'зелёный' (zelyonyy) sounds very much like the English surname 'Greene'.")

    print("\nStep 4: We need to find an English poet with this surname. Robert Greene (1558–1592) was a prominent English poet and playwright.")

    final_surname = "Greene"
    print("\nTherefore, the surname of the English poet that appears in the description is Greene.")

# The following lines simulate the final equation by showing the numbered steps that lead to the answer.
print("Final logical chain:")
print("1 (Vienna's Green Boulevards) leads to")
print("2 (Russian description 'зелёный') which sounds like")
print("3 (The surname of English poet Robert Greene).")

# Using sys.stdout.write to prevent an extra newline and have the answer appear on a fresh line
# as if it's the only output, per some platform requirements.
# sys.stdout.write(final_surname)

solve_riddle()