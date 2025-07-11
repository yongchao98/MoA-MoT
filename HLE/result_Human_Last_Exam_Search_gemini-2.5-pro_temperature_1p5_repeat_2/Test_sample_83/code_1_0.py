def solve_multiple_knapsack_greedy():
    """
    Solves the multiple knapsack problem using a greedy algorithm based on value/weight ratio.
    This is a heuristic and may not produce the optimal solution.
    """
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    knapsack_capacities = [40, 120, 200]

    num_items = len(values)
    num_knapsacks = len(knapsack_capacities)

    # Create a list of items with their value, weight, ratio, and original index
    items = []
    for i in range(num_items):
        # To avoid division by zero if weight is 0
        ratio = values[i] / weights[i] if weights[i] != 0 else values[i]
        items.append({'value': values[i], 'weight': weights[i], 'ratio': ratio, 'id': i})

    # Sort items by value-to-weight ratio in descending order
    items.sort(key=lambda x: x['ratio'], reverse=True)

    # Track knapsack loads and contents
    knapsack_loads = [0] * num_knapsacks
    knapsack_contents = [[] for _ in range(num_knapsacks)]
    total_value = 0
    
    # Greedily fill knapsacks
    for item in items:
        for i in range(num_knapsacks):
            # Check if the item fits in the current knapsack
            if knapsack_loads[i] + item['weight'] <= knapsack_capacities[i]:
                # Add item to the knapsack
                knapsack_loads[i] += item['weight']
                knapsack_contents[i].append(item['id'])
                total_value += item['value']
                # Move to the next item
                break
                
    return total_value

# Since this is a direct calculation, I will just call the function.
# The result will be a numerical value.
final_value = solve_multiple_knapsack_greedy()
print(final_value)
