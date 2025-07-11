import sys

def solve_treewidth_bound():
    """
    This function takes treewidths of two graphs H and G, and the size
    of their vertex intersection, then calculates the tight upper bound on
    the treewidth of their union F.

    The user is expected to provide three integer values as command-line arguments:
    t_H: treewidth of graph H
    t_G: treewidth of graph G
    k: number of vertices in the intersection V(H) intersection V(G)
    """
    # Prompt the user if inputs are not provided via command line
    if len(sys.argv) < 4:
        try:
            t_H_str = input("Enter the treewidth of graph H (t_H): ")
            t_G_str = input("Enter the treewidth of graph G (t_G): ")
            k_str = input("Enter the number of intersection vertices (k): ")
            t_H = int(t_H_str)
            t_G = int(t_G_str)
            k = int(k_str)
        except (ValueError, TypeError):
            print("Invalid input. Please enter integer values.")
            return
    else:
        try:
            t_H = int(sys.argv[1])
            t_G = int(sys.argv[2])
            k = int(sys.argv[3])
        except (ValueError, IndexError):
            print("Usage: python your_script.py <t_H> <t_G> <k>")
            print("Please provide three integer arguments.")
            return
    
    # Calculate the maximum of the two treewidths
    max_t = max(t_H, t_G)
    
    # Calculate the upper bound on the treewidth of F
    bound = max_t + k
    
    # Print the equation step by step
    print(f"Given values: t_H = {t_H}, t_G = {t_G}, k = {k}")
    print(f"The tight upper bound for the treewidth of F is given by the formula: max(t_H, t_G) + k")
    print(f"Calculation: max({t_H}, {t_G}) + {k} = {max_t} + {k} = {bound}")
    
if __name__ == '__main__':
    solve_treewidth_bound()
