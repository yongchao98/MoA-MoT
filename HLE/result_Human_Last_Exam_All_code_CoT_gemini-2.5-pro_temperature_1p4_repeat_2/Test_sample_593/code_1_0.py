def solve_treewidth_bound():
    """
    Calculates and prints a tight upper bound on the treewidth of the union of two graphs.
    
    The user is prompted to enter the treewidth of graph H (t_H),
    the treewidth of graph G (t_G), and the number of shared vertices (k).
    """
    try:
        t_H_str = input("Enter the treewidth of H (t_H): ")
        t_H = int(t_H_str)

        t_G_str = input("Enter the treewidth of G (t_G): ")
        t_G = int(t_G_str)

        k_str = input("Enter the number of shared vertices (k): ")
        k = int(k_str)

        if k <= 0:
            print("The number of shared vertices (k) must be positive.")
            return

        # The tight upper bound for the treewidth of the combined graph F is max(t_H, t_G, k-1).
        k_minus_1 = k - 1
        result = max(t_H, t_G, k_minus_1)
        
        # Print the equation with all the numbers, as requested.
        print(f"The tight upper bound for the treewidth of F is calculated as:")
        print(f"max(t_H, t_G, k-1)")
        print(f"= max({t_H}, {t_G}, {k}-1)")
        print(f"= max({t_H}, {t_G}, {k_minus_1})")
        print(f"= {result}")

    except ValueError:
        print("Invalid input. Please enter integer values.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_treewidth_bound()