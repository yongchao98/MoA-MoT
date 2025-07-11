import sys

def calculate_max_rank():
    """
    Calculates and explains the maximal rank of the complementary channel
    of a Pauli channel on a d-dimensional system.
    """
    try:
        # Prompt the user for the dimension 'd' of the qudit.
        d_input = input("Enter the dimension of the qudit (d >= 2): ")
        d = int(d_input)

        if d < 2:
            print("Error: The dimension d must be an integer greater than or equal to 2.", file=sys.stderr)
            return

        # The maximal rank is d^2.
        max_rank = d ** 2

        # Print the explanation and the result.
        print("\n--- Derivation of Maximal Rank ---")
        print(f"For a Pauli channel on a d={d} dimensional system (a qudit):")
        print("\nStep 1: The channel is defined by a sum over the d^2 generalized Pauli operators.")
        print(f"   There are d*d = {d}*{d} = {d**2} such operators.")
        
        print("\nStep 2: The rank of the complementary channel equals the rank of the Gram matrix G of the channel's Kraus operators.")
        print("   The Kraus operators are proportional to the Pauli operators.")
        
        print("\nStep 3: Due to the orthogonality of Pauli operators, the Gram matrix G is diagonal.")
        print("   Each diagonal element is proportional to one of the probabilities defining the channel.")

        print("\nStep 4: The rank of this diagonal matrix is the number of its non-zero diagonal entries.")
        print("   This is equal to the number of non-zero probabilities used to define the channel.")

        print("\nStep 5: To maximize the rank, we must maximize the number of non-zero probabilities.")
        print(f"   We can choose all {d**2} probabilities to be non-zero (e.g., p_k = 1/{d**2} for all k).")
        
        print("\n--- Conclusion ---")
        print(f"The maximal rank is achieved when all {d**2} Kraus operators are used.")
        print(f"Therefore, the maximal rank is {d} ^ 2 = {max_rank}.")

    except ValueError:
        print("Error: Invalid input. Please enter an integer.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    calculate_max_rank()