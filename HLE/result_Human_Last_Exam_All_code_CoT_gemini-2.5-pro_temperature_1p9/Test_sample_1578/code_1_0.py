import sys

def solve_mirror_puzzle():
    """
    This script simulates a test to determine if a mirror's reflection
    is real or an illusion created by an Oni (demon).

    The test is based on a core principle of optics: a reflection requires
    a clear line of sight. By placing an opaque barrier between yourself
    and the mirror, you can verify this principle.

    - A real mirror will no longer show your reflection.
    - An Oni creating an illusion might not account for the barrier and
      continue to project your image, revealing the deception.
    """

    # --- Scenario Setup ---
    # We can model the scene with simple 1D positions.
    your_position = 0
    mirror_position = 10
    barrier_position = 5  # Placed between you and the mirror

    print("The Japanese Mirror and Oni Puzzle\n")
    print("Testing the principle with a known REAL mirror...")
    print("-----------------------------------------------")
    print(f"Your Position: {your_position}")
    print(f"Mirror Position: {mirror_position}")
    print(f"Opaque Barrier Position: {barrier_position}")

    # For a real mirror, the line of sight is blocked.
    is_line_of_sight_blocked = (your_position < barrier_position < mirror_position)

    if is_line_of_sight_blocked:
        print("\nResult: Your reflection is NOT visible.")
        print("This is the expected behavior for a real mirror.")
    else:
        # This part of the conditional would run if the barrier wasn't between you and the mirror.
        print("\nResult: Your reflection is visible.")


    print("\nNow, testing the SUSPICIOUS mirror...")
    print("-------------------------------------")
    print(f"You are at position: {your_position}")
    print(f"The suspicious mirror is at position: {mirror_position}")
    print(f"You place the opaque barrier at position: {barrier_position}")

    # We assume the Oni does not respect the physical line-of-sight barrier.
    # It sees you and projects an illusion of you "in" the mirror.
    oni_shows_illusion = True

    if is_line_of_sight_blocked and oni_shows_illusion:
        print("\nResult: Your image IS STILL VISIBLE in the mirror!")
        print("\nConclusion: This is not a reflection! The physical line of sight is blocked.")
        print("The figure in the mirror is an illusion created by an Oni.")
    else:
        print("\nResult: The test is inconclusive or the Oni is very clever.")
        print("However, in this scenario, the Oni has revealed itself.")

if __name__ == '__main__':
    solve_mirror_puzzle()
