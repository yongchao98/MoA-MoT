# This script determines and prints the significant contributors to the observation of
# dual light/alpha chain expression in B and T cells from scRNA-seq data.

# Based on the analysis, the significant contributors are identified.
# For B cells, these include technical artifacts (doublets, ambient RNA) and
# biological mechanisms (functional allelic inclusion and receptor editing).
b_cell_contributors = [1, 2, 3, 5]

# For T cells, these include technical artifacts and the common biological
# phenomenon of dual alpha chain expression. Autoreactivity is not a major
# factor due to negative selection.
t_cell_contributors = [1, 2, 3]

# The lists are formatted into strings with numbers separated by commas.
# The str() function is used to convert each number to a string.
# The ",".join() method concatenates them with a comma in between.
b_cell_str = ",".join(map(str, b_cell_contributors))
t_cell_str = ",".join(map(str, t_cell_contributors))

# The final output string is constructed according to the specified format:
# (b_cell_list), (t_cell_list)
final_answer_string = f"({b_cell_str}), ({t_cell_str})"

# The final formatted string is printed to the console.
print(final_answer_string)