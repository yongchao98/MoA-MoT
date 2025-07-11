import sys
import time

def simulate_cosmo_loco_output():
    """
    This function simulates the output of the given Scheme code.
    The Scheme code prints 1, then 1 again, and then enters an
    infinite loop printing 2.
    """
    try:
        # The first `(write 1)`
        print(1, end=' ', flush=True)

        # The second `(write 1)`
        print(1, end=' ', flush=True)

        # The infinite loop printing `(write 2)`
        while True:
            print(2, end=' ', flush=True)
            # We add a small delay to prevent the terminal from being overwhelmed
            # and to make the infinite stream of 2s observable.
            time.sleep(0.5)
            
    except KeyboardInterrupt:
        # Allow the user to stop the infinite loop with Ctrl+C
        print("\nSimulation stopped by user.")

if __name__ == "__main__":
    simulate_cosmo_loco_output()