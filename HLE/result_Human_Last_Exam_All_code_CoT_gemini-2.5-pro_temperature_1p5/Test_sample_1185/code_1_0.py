def main():
    """
    This script enumerates and describes the different types of stable reduction
    for curves of genus 2, and prints the total count.
    """

    # List of known types of stable reduction for curves of genus 2
    # Each entry contains a description and parameters for the genus formula:
    # p_a = (sum_of_geometric_genera - num_components + 1) + num_nodes
    stable_reductions = [
        {
            "description": "An irreducible curve with one node",
            "geo_genera": [1],
            "num_components": 1,
            "num_nodes": 1
        },
        {
            "description": "An irreducible curve with two nodes",
            "geo_genera": [0],
            "num_components": 1,
            "num_nodes": 2
        },
        {
            "description": "Two elliptic curves meeting at one point",
            "geo_genera": [1, 1],
            "num_components": 2,
            "num_nodes": 1
        },
        {
            "description": "An elliptic curve and a rational curve with a node, meeting at another point",
            "geo_genera": [1, 0],
            "num_components": 2,
            "num_nodes": 2
        },
        {
            "description": "Two rational curves meeting at three points",
            "geo_genera": [0, 0],
            "num_components": 2,
            "num_nodes": 3
        }
    ]

    print("The different types of stable reduction for a curve of genus 2 are listed below.\n")
    print("For each type, we verify that the arithmetic genus is 2 using the formula:")
    print("p_a = (sum of geometric genera - number of components + 1) + number of nodes\n")
    print("-" * 50)

    for i, reduction_type in enumerate(stable_reductions):
        description = reduction_type["description"]
        geo_genera = reduction_type["geo_genera"]
        k = reduction_type["num_components"]
        delta = reduction_type["num_nodes"]
        
        sum_g = sum(geo_genera)
        # p_a = (sum_g - k + 1) + delta
        p_a = (sum_g - k + 1) + delta

        print(f"Type {i+1}: {description}")
        print(f"  - Number of components (k): {k}")
        print(f"  - Geometric genera of components (g_i): {geo_genera}")
        print(f"  - Number of nodes (delta): {delta}")
        print(f"  - Arithmetic Genus Equation: p_a = (sum({geo_genera}) - {k} + 1) + {delta}")
        print(f"  - Calculation: p_a = ({sum_g} - {k} + 1) + {delta} = {p_a}")
        print("-" * 50)

    count = len(stable_reductions)
    print(f"\nBased on this classification, there are {count} different combinatorial types of stable reduction for curves of genus 2.")

if __name__ == "__main__":
    main()
<<<5>>>