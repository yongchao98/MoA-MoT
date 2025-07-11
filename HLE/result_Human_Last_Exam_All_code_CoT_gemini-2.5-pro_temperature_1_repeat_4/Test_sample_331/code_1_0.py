def solve_chow_group_problem():
    """
    This function prints the computed pairs (m(X), M(X)) for the four given varieties.
    The values are pre-computed based on principles of algebraic geometry.

    - X_1 (genus 2 curve): m(X_1)=2, M(X_1)=3
    - X_2 (general genus 7 curve): m(X_2)=8, M(X_2)=8
    - X_3 (Enriques surface): m(X_3)=2, M(X_3)=2
    - X_4 (Grassmannian G(3,6)): m(X_4)=1, M(X_4)=1
    """

    # The computed pairs (m(X), M(X))
    result_x1 = (2, 3)
    result_x2 = (8, 8)
    result_x3 = (2, 2)
    result_x4 = (1, 1)

    # Format the output string as requested
    # The output needs to show each number in the final result.
    output_string = f"({result_x1[0]}, {result_x1[1]}), ({result_x2[0]}, {result_x2[1]}), ({result_x3[0]}, {result_x3[1]}), ({result_x4[0]}, {result_x4[1]})"
    
    print(output_string)

solve_chow_group_problem()