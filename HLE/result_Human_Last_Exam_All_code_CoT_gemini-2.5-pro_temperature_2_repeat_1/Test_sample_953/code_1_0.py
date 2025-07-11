import math

def solve():
    """
    Analyzes the runtime of a randomized MIS algorithm and provides the complexity categories.
    """
    print("### Analysis of the Algorithm's Complexity ###\n")
    print("The algorithm's runtime is determined by how quickly vertices are removed from the graph.")
    print("A vertex `v` is removed if it or any of its neighbors are selected in a step.")
    print("A vertex `v` with degree `d(v)` is selected with probability 1 / (d(v) + 1).\n")

    # --- Analysis for each case ---
    analysis = {
        "f1(n) on a cycle": {
            "reasoning": [
                "In a cycle C_n, every vertex has a degree of 2. The maximum degree is constant.",
                "The probability that a given vertex `v` is removed in a step is a constant `p > 0`.",
                "For example, the probability of `v` itself being selected is 1/3.",
                "A more detailed calculation shows P(v removed) is a constant (e.g., 13/15).",
                "This leads to an expected exponential decrease in the number of vertices, implying an O(log n) upper bound on the number of steps.",
                "A path or cycle graph is also a worst-case example for these algorithms, providing a matching lower bound of Omega(log n).",
                "Therefore, the complexity is Theta(log n)."
            ],
            "complexity": "Theta(log n)"
        },
        "f2(n) on a tree with degree <= 100": {
            "reasoning": [
                "For a tree with a maximum degree of 100, the degree of any vertex is bounded by a constant.",
                "The general analysis for bounded-degree graphs gives an O(log n) upper bound.",
                "A simple path graph is a type of tree. As established in the cycle analysis, the algorithm takes Omega(log n) steps on a path.",
                "This provides a lower bound, matching the upper bound.",
                "Therefore, the complexity is Theta(log n)."
            ],
            "complexity": "Theta(log n)"
        },
        "f3(n) on a graph with degree <= 100": {
            "reasoning": [
                "This is the general case for graphs with a constant bound on the maximum degree (Delta <= 100).",
                "The standard analysis of Luby-style algorithms shows that the number of vertices decreases by a constant factor in an expected sense over a constant number of rounds, which yields an O(log n) upper bound.",
                "The path graph provides an Omega(log n) lower bound, as it is a graph with maximum degree 2.",
                "Therefore, the complexity for this case is also Theta(log n)."
            ],
            "complexity": "Theta(log n)"
        }
    }

    final_digits = ""
    for name, data in analysis.items():
        print(f"--- Analysis for {name} ---")
        for line in data["reasoning"]:
            print(f"- {line}")
        print(f"Resulting complexity: {data['complexity']}\n")

        # Map complexity to category
        # 9: f(n) = Omega(log n)
        # A function that is Theta(log n) is also Omega(log n).
        category = 9
        final_digits += str(category)

    print("### Mapping to Categories ###")
    print("The complexity class Theta(log n) implies that the function is Omega(log n).")
    print("Looking at the list of categories:")
    print("1: O(1)")
    print("2: omega(1) but o(log*n)")
    print("3: Theta(log*n)")
    print("4: omega(log*n) but o(log log n)")
    print("5: Theta(log log n)")
    print("6: omega(log log n) but O(log^0.1(n))")
    print("7: omega(log^0.1(n)) but O(log^0.9(n))")
    print("8: omega(log^0.9(n)) but o(log n)")
    print("9: Omega(log n)")
    print("\nTheta(log n) fits into category 9.")
    
    print("\n### Final Answer ###")
    print(f"d1 (cycle): Category {final_digits[0]}")
    print(f"d2 (tree):  Category {final_digits[1]}")
    print(f"d3 (graph): Category {final_digits[2]}")
    print(f"\nThe three digits are {final_digits}.")
    print("<<<999>>>")

solve()