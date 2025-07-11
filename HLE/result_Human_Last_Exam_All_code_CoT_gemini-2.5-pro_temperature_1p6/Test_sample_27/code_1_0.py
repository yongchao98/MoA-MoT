def solve_bcr_tcr_question():
    """
    This function determines the significant contributors to observing
    dual light/alpha chains in single-cell B/T cell data and prints the answer.
    """

    # For B cells, technical artifacts are significant. The primary biological
    # mechanism is receptor editing to resolve autoreactivity (5), which implies
    # the original chain's protein is non-functional or not on the surface (4, 6).
    # True dual-functional expression (3) also occurs at significant levels.
    # Therefore, all listed mechanisms are considered significant.
    b_cell_causes = [1, 2, 3, 4, 5, 6]

    # For T cells, technical artifacts are also significant. The primary biological
    # cause is the lack of allelic exclusion at the alpha locus (3), which is very
    # common. Consequently, one chain might be less functional or fail to reach
    # the surface (4, 6). However, cells with an autoreactive chain (5) are
    # strongly selected against in the thymus, making this an insignificant
    # contributor to the observed naive T cell population.
    t_cell_causes = [1, 2, 3, 4, 6]

    # Format the cause lists into the required string format "(n,n,...)"
    # The code explicitly constructs the string from each number in the lists.
    b_cell_str_list = []
    for num in b_cell_causes:
        b_cell_str_list.append(str(num))
    b_cell_str = "(" + ",".join(b_cell_str_list) + ")"

    t_cell_str_list = []
    for num in t_cell_causes:
        t_cell_str_list.append(str(num))
    t_cell_str = "(" + ",".join(t_cell_str_list) + ")"

    # Combine the two parts into the final answer string
    final_answer = b_cell_str + ", " + t_cell_str
    
    print(final_answer)

solve_bcr_tcr_question()
print("<<<" + "(1,2,3,4,5,6), (1,2,3,4,6)" + ">>>")