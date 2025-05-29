def greedy_knapsack(profit, weight, max_weight):
    if max_weight <= 0:
        raise ValueError("max_weight must be greater than zero.")
    if len(profit) != len(weight):
        raise IndexError("The length of profit and weight must be the same.")
    if any(p < 0 for p in profit):
        raise ValueError("Profit cannot be negative.")
    if any(w < 0 for w in weight):
        raise ValueError("Weight cannot be negative.")
    
    # Greedy knapsack algorithm
    items = list(zip(profit, weight))
    items.sort(key=lambda x: x[0] / x[1], reverse=True)
    
    total_profit = 0
    remaining_weight = max_weight
    
    for p, w in items:
        if remaining_weight >= w:
            total_profit += p
            remaining_weight -= w
        else:
            total_profit += (p / w) * remaining_weight
            break
    
    return {"total_profit": total_profit}

# Given input
input_data = {'profit': [12, 80, 86, 96, 82, 66, 5, 81, 48], 'weight': [57, 68, 14, 88, 59, 77, 37, 49, 95], 'max_weight': 158}

# Calculate the result
result = greedy_knapsack(input_data['profit'], input_data['weight'], input_data['max_weight'])
print(result)