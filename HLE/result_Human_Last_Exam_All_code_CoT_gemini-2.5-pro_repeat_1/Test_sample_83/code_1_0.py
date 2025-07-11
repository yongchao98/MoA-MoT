import sys

class MultipleKnapsackSolver:
    """
    Solves the Multiple Knapsack Problem where each item can be used at most once,
    using dynamic programming with memoization.
    """

    def __init__(self, values, weights, capacities):
        self.values = values
        self.weights = weights
        self.capacities = capacities
        self.num_items = len(values)
        self.num_knapsacks = len(capacities)
        self.memo = {}
        # The recursion depth will be at most the number of items.
        # Python's default limit is usually high enough, but we can increase it if needed.
        required_recursion_limit = self.num_items + 5
        if sys.getrecursionlimit() < required_recursion_limit:
            sys.setrecursionlimit(required_recursion_limit)

    def _solve_recursive(self, item_index, current_capacities):
        """
        Recursively finds the maximum value using memoization.

        Args:
            item_index (int): Index of the current item being considered.
            current_capacities (tuple): A tuple of remaining capacities for the knapsacks.

        Returns:
            int: The maximum value achievable from this state.
        """
        # Base case: If we have considered all items, no more value can be added.
        if item_index == self.num_items:
            return 0

        # Use the state (item index, capacities tuple) as the key for memoization.
        state = (item_index, current_capacities)
        if state in self.memo:
            return self.memo[state]

        # --- Recursive Step ---

        # Option 1: Don't place the current item in any knapsack.
        # The value is what we can get from the remaining items.
        max_val = self._solve_recursive(item_index + 1, current_capacities)

        # Option 2: Try to place the current item in each knapsack.
        current_value = self.values[item_index]
        current_weight = self.weights[item_index]

        for k in range(self.num_knapsacks):
            # If the item fits in the current knapsack
            if current_weight <= current_capacities[k]:
                # Create a new capacities tuple for the recursive call.
                new_caps_list = list(current_capacities)
                new_caps_list[k] -= current_weight
                new_capacities_tuple = tuple(new_caps_list)
                
                # Calculate the value for this choice.
                val_if_placed = current_value + self._solve_recursive(item_index + 1, new_capacities_tuple)
                
                # Update max_val if this choice is better.
                max_val = max(max_val, val_if_placed)

        # Store the result in the memoization table and return it.
        self.memo[state] = max_val
        return max_val

    def _reconstruct_solution(self):
        """
        Backtracks through the memoization table to find the items placed in each knapsack.
        
        Returns:
            tuple: Contains two lists of lists (knapsack_values, knapsack_weights)
                   detailing the items assigned to each knapsack.
        """
        knapsacks_values = [[] for _ in range(self.num_knapsacks)]
        knapsacks_weights = [[] for _ in range(self.num_knapsacks)]
        
        current_capacities = tuple(self.capacities)
        
        # Iterate through each item to determine where it was placed.
        for i in range(self.num_items):
            state = (i, current_capacities)
            
            # The value if we skip item i is found in the memoized result for the next item.
            val_if_skipped = self.memo.get((i + 1, current_capacities), 0)
            
            # If the optimal value for the current state is the same as if we skipped the item,
            # then the item was not taken.
            if self.memo[state] == val_if_skipped:
                continue
            
            # Otherwise, the item was taken. We need to find which knapsack it was placed in.
            current_value = self.values[i]
            current_weight = self.weights[i]
            
            for k in range(self.num_knapsacks):
                if current_weight <= current_capacities[k]:
                    new_caps_list = list(current_capacities)
                    new_caps_list[k] -= current_weight
                    new_capacities_tuple = tuple(new_caps_list)

                    # Check if placing item i in knapsack k leads to the optimal value.
                    val_if_placed = current_value + self.memo.get((i + 1, new_capacities_tuple), 0)
                    
                    if self.memo[state] == val_if_placed:
                        knapsacks_values[k].append(current_value)
                        knapsacks_weights[k].append(current_weight)
                        
                        # Update capacities for the next iteration.
                        current_capacities = new_capacities_tuple
                        break  # Move to the next item once placement is found.

        return knapsacks_values, knapsacks_weights
        
    def solve_and_print(self):
        """
        Executes the solver and prints the formatted solution.
        """
        # Step 1: Calculate the maximum value by populating the memoization table.
        initial_capacities = tuple(self.capacities)
        max_total_value = self._solve_recursive(0, initial_capacities)

        # Step 2: Reconstruct the solution to find item assignments.
        knapsacks_values, knapsacks_weights = self._reconstruct_solution()
        
        # Step 3: Print the results.
        print(f"Maximum total value: {max_total_value}\n")
        
        all_item_values = []
        for k in range(self.num_knapsacks):
            total_k_value = sum(knapsacks_values[k])
            total_k_weight = sum(knapsacks_weights[k])
            all_item_values.extend(knapsacks_values[k])
            
            print(f"Knapsack {k+1} (Capacity: {self.capacities[k]})")
            print(f"  Total Weight: {total_k_weight}")
            print(f"  Total Value: {total_k_value}")
            
            if knapsacks_values[k]:
                value_str = " + ".join(map(str, sorted(knapsacks_values[k], reverse=True)))
                print(f"  Values of items: {value_str}\n")
            else:
                print(f"  Values of items: (empty)\n")

        print("-" * 30)
        
        # Print the final equation for the total value.
        print("Final equation for the total value:")
        if all_item_values:
            final_equation_str = " + ".join(map(str, sorted(all_item_values, reverse=True)))
            print(f"{final_equation_str} = {max_total_value}")
        else:
            print("0 = 0")


# --- Problem Data ---
item_values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
item_weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
knapsack_capacities = [40, 120, 200]

# --- Execution ---
solver = MultipleKnapsackSolver(item_values, item_weights, knapsack_capacities)
solver.solve_and_print()
<<<694>>>