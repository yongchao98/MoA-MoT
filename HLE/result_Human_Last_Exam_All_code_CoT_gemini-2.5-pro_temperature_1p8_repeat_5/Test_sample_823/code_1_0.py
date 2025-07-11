def explain_reasoning():
    """
    This script explains the reasoning for choosing option E.
    The argument rests on the Excluded Grid Minor Theorem, a cornerstone of graph theory.
    """
    
    reasoning_steps = [
        "1. Definition of Unbounded Treewidth:",
        "   The class C has 'unbounded treewidth'. This means for any integer W, there is a graph G in C with treewidth(G) > W.",
        "",
        "2. The Excluded Grid Minor Theorem:",
        "   A central result by Robertson and Seymour states that large treewidth implies a large grid minor. Formally, for any k, there exists a treewidth W(k) such that any graph with treewidth >= W(k) contains a k-by-k grid as a minor. This remains true for graphs with a bounded maximum degree.",
        "",
        "3. From Grid Minor to Clique Minor:",
        "   A k-by-k grid graph is known to contain a clique of size proportional to k as a minor (e.g., K_floor(k/2)). A clique is a graph where every vertex is connected to every other vertex. A graph H is a minor of G if H can be obtained from a subgraph of G by contracting edges.",
        "",
        "4. Conclusion:",
        "   - Since C has unbounded treewidth, it contains graphs with treewidth large enough to guarantee any k-by-k grid minor.",
        "   - The presence of arbitrarily large grid minors implies the presence of arbitrarily large clique minors.",
        "   - Therefore, the graphs in C must contain clique-minors of unbounded size.",
        "   This matches statement E."
    ]
    
    print("Reasoning for selecting option E:")
    for step in reasoning_steps:
        print(step)

explain_reasoning()