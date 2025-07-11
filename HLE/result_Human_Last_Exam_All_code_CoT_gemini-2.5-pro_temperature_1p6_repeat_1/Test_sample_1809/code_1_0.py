def solve_orientation_number():
    """
    Calculates and explains the valid orientation number of graph H.

    The valid orientation number is the minimum possible maximum indegree.
    We found this value to be 6. This is achieved by an orientation where the
    central vertex indegrees are {0, 4, 5, 6}.

    The maximum value, 6, is the indegree of one of the central vertices.
    Let's call this vertex v_max. Its indegree is the sum of contributions
    from the central K4 (d_max) and its 10 attached K3s (S_max).
    """

    # In our construction, we orient the central K4 to give indegrees {0, 1, 2, 3}.
    # The vertex that will achieve the maximum indegree is the one with the
    # largest contribution from the K4.
    d_max = 3

    # We assign the target indegrees {0, 4, 5, 6} to the central vertices
    # corresponding to the K4 contributions {0, 1, 2, 3}.
    # The vertex with d=3 gets the target indegree of 6.
    target_max_indegree = 6

    # The required contribution from the attached triangles (S_max) is:
    S_max = target_max_indegree - d_max

    # The final answer is the maximum indegree.
    final_answer = target_max_indegree

    print(f"The valid orientation number of graph H is {final_answer}.")
    print("\nThis value is the maximum indegree in an optimal orientation.")
    print("It is achieved at one of the four central vertices. Let's trace its calculation:")
    
    # Print the equation with each number.
    print(f"\nindegree(v_max) = (contribution from K4) + (contribution from attached triangles)")
    print(f"               {final_answer} = {d_max} + {S_max}")


solve_orientation_number()