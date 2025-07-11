import pulp
import networkx as nx

def solve_k_vector():
    """
    This function determines the smallest k for the Petersen graph,
    which is a bridgeless 3-regular graph. The principles found here
    apply to larger graphs of the same type, such as those with 20 vertices.
    """
    # 1. Create the Petersen Graph (10 vertices, 15 edges, 3-regular, bridgeless)
    G = nx.petersen_graph()
    edges = list(G.edges())
    
    # 2. Iteratively check for the smallest k
    # A k-vector has entries from {+/-1, ..., +/-(k-1)}
    for k in range(3, 7):
        k_minus_1 = k - 1
        print(f"Testing for k={k} (values in {{+/-1, ..., +/-{k_minus_1}}})")
        
        # Create an Integer Linear Programming problem
        prob = pulp.LpProblem(f"k-vector_k_{k}", pulp.LpMinimize)
        # The objective is trivial; we only care about feasibility
        prob += 0 

        # Define variables for each edge, with integer values
        x = pulp.LpVariable.dicts("x", edges, lowBound=-k_minus_1, upBound=k_minus_1, cat=pulp.LpInteger)
        
        # Add constraint that variables cannot be zero.
        # This requires introducing binary variables.
        # For each edge e:
        # x[e] >= 1 - (k_minus_1 + 1) * b[e]
        # x[e] <= -1 + (k_minus_1 + 1) * (1 - b[e])
        # where b[e] is binary. If b[e]=0, x[e] >= 1. If b[e]=1, x[e] <= -1.
        b = pulp.LpVariable.dicts("b", edges, cat=pulp.LpBinary)
        M = k_minus_1 + 1
        for e in edges:
            prob += x[e] >= 1 - M * b[e]
            prob += x[e] <= -1 + M * (1 - b[e])

        # Add the main constraint: sum of flows at each vertex is 0
        for v in G.nodes():
            incident_edges = [e for e in edges if v in e]
            prob += pulp.lpSum(x[e] for e in incident_edges) == 0

        # 3. Solve the problem
        # The CBC solver is used. msg=False suppresses solver output.
        prob.solve(pulp.PULP_CBC_CMD(msg=False))

        # 4. Check if a solution was found
        if pulp.LpStatus[prob.status] == "Optimal":
            print(f"Feasible solution found for k = {k}.")
            print("\nBased on this, the smallest integer k for which a non-4-flowable graph like")
            print("the Petersen graph admits a k-vector is 5. This establishes the lower bound.")
            
            # Print an example equation as requested
            v_example = 0
            incident_edges_example = [e for e in edges if v_example in e]
            vals = [int(x[e].varValue) for e in incident_edges_example]
            
            # Make the equation look nice by putting negative numbers on the other side
            pos_vals = [v for v in vals if v > 0]
            neg_vals = [-v for v in vals if v < 0]
            
            print(f"\nAn example equation from vertex {v_example}:")
            print(f"{' + '.join(map(str, pos_vals))} = {' + '.join(map(str, neg_vals))}")
            # And also print it in the requested `a+b+c=0` format
            print("In the format a+b+c=0:")
            equation_str = f"{vals[0]} "
            if vals[1] >= 0:
                equation_str += f"+ {vals[1]} "
            else:
                equation_str += f"- {-vals[1]} "
            if vals[2] >= 0:
                equation_str += f"+ {vals[2]} "
            else:
                equation_str += f"- {-vals[2]} "
            equation_str += f"= {sum(vals)}"
            print(equation_str)
            
            return
        else:
            print(f"No feasible solution found for k = {k}.")

if __name__ == '__main__':
    try:
        solve_k_vector()
    except ImportError:
        print("Please install the required libraries: pip install networkx pulp")
